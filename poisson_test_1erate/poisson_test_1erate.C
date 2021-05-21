#include "TPaveLabel.h"
#include "TLine.h"
#include <random>
#include <ctime>
#include <cstdlib>
#include "fitsio.h"
#include <stdio.h>
#include <stdlib.h>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <iostream>
#include <chrono>
#include <sstream>
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TStopwatch.h"
#include <TROOT.h>
#include <TStyle.h>
#include "TF1.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TGraphErrors.h"
#include "TLeafC.h"
#include "TStyle.h"
#include "TChain.h"
#include "assert.h"

using namespace std;

#define ARRSIZE(arr) (sizeof(arr)/sizeof(*(arr)))

// Masfile configuration variables // 
int NROW;int NCOL;int NBINROW;int NBINCOL;int CCDNCOL;int CCDNROW;int CCDNPRES;int RUNID;int LTANAME;

// compare protocol
int compare(const void *a, const void *b) {
    double x1 = *(const double*)a;
    double x2 = *(const double*)b;
    if (x1 > x2) return  1;
    if (x1 < x2) return -1;
    // x1 and x2 are equal; compare y's
    double y1 = *(((const double*)a)+1);
    double y2 = *(((const double*)b)+1);
    if (y1 > y2) return  1;
    if (y1 < y2) return -1;
    return 0;
}

void concatenate_doubles(double result[][3],double a[][3], double b[][3], int size_a, int size_b)
{
    for (size_t i = 0; i < size_a; i++)
    {
        result[i][0]=a[i][0];
        result[i][1]=a[i][1];
        result[i][2]=a[i][2];
    }
    for (size_t i = 0; i < size_b; i++)
    {
        result[i+size_a][0]=b[i][0];
        result[i+size_a][1]=b[i][1];
        result[i+size_a][2]=b[i][2];
    }
    // qsort(result, ARRSIZE(result), sizeof(*result), compare);
    

    return;
}

void compute_largest_diff(double ladder[][3],vector <double> pvalues, double &largest_diff)
{
    double diff,last_ON,last_OFF=0;
    for (Int_t i=0;i<100+pvalues.size();i++)
    {
        if (ladder[i][2]==0) last_OFF = ladder[i][1];
        if (ladder[i][2]==1) last_ON = ladder[i][1];
        diff=abs(last_ON-last_OFF);
        if (largest_diff<diff) largest_diff=diff;
    }
    return;
}



// weighted_average
double weighted_average( const vector <float> std, const vector <float> scores, std::size_t n )
{
    // see: https://en.wikipedia.org/wiki/Weighted_arithmetic_mean#Mathematical_definition

    // compute the sum of weights
    float sum_of_weights = 0 ;
    for( std::size_t i = 0 ; i < n ; ++i ) sum_of_weights += 1/pow(std[i],2) ;

    // compute the total of each score multiplied by its weight
    float weighted_total = 0 ;
    for( std::size_t i = 0 ; i < n ; ++i ) weighted_total += scores[i] * (1/pow(std[i],2)) ;

    // divide the weighted_total by the sum_of_weights to get the weighted average
    return weighted_total / sum_of_weights ;
}

// Standard deviation //
double stddev(std::vector<float> const & func)
{
    double mean = std::accumulate(func.begin(), func.end(), 0.0) / func.size();
    double sq_sum = std::inner_product(func.begin(), func.end(), func.begin(), 0.0,
        [](double const & x, double const & y) { return x + y; },
        [mean](double const & x, double const & y) { return (x - mean)*(y - mean); });
    return std::sqrt(sq_sum / func.size());
}

//////////////////Read maskfile configuration ////////////////////////////

void readConfig(TFile *file)
{
	TTree* config =nullptr; file->GetObject("headerTree_0",config); config->GetEvent(0);
	char* NROWchar;NROWchar=((TLeafC *) config->GetBranch("NROW")->GetLeaf("string"))->GetValueString();istringstream buffer(NROWchar); buffer >> NROW;
    cout << "NROW: " << NROW << endl; 
	char* NCOLchar;NCOLchar=((TLeafC *) config->GetBranch("NCOL")->GetLeaf("string"))->GetValueString();istringstream buffer2(NCOLchar); buffer2 >> NCOL; 
    cout << "NCOL: " << NCOL << endl;
	char* CCDNROWchar;CCDNROWchar=((TLeafC *) config->GetBranch("CCDNROW")->GetLeaf("string"))->GetValueString();istringstream buffer3(CCDNROWchar); buffer3 >> CCDNROW; 
    cout << "CCDNROW: " << CCDNROW << endl;
	char* CCDNCOLchar;CCDNCOLchar=((TLeafC *) config->GetBranch("CCDNCOL")->GetLeaf("string"))->GetValueString();istringstream buffer4(CCDNCOLchar); buffer4 >> CCDNCOL;
    cout << "CCDNCOL: " << CCDNCOL << endl;
    CCDNPRES=7;
    NBINROW=1;
    NBINCOL=1;
    RUNID=1;
    LTANAME=1;
	// try {char CCDNPRESchar;CCDNPRESchar=((TLeafC *) config->GetBranch("CCDNPRES")->GetLeaf("string"))->GetValueString();istringstream buffer5(CCDNPRESchar); buffer5 >> CCDNPRES; throw 501;} catch (...){CCDNPRES=7;}
    // try {char* NBINROWchar;NBINROWchar=((TLeafC *) config->GetBranch("NBINROW")->GetLeaf("string"))->GetValueString();istringstream buffer6(NBINROWchar); buffer6 >> NBINROW; 
    // throw 502;} catch (...){NBINROW=1;}
	// try {char* NBINCOLchar;NBINCOLchar=((TLeafC *) config->GetBranch("NBINCOL")->GetLeaf("string"))->GetValueString();istringstream buffer7(NBINCOLchar); buffer7 >> NBINCOL; throw 503;} catch (...){NBINCOL=1;}
	// try{char* RUNIDchar;RUNIDchar=((TLeafC *) config->GetBranch("RUNID")->GetLeaf("string"))->GetValueString();istringstream buffer8(RUNIDchar); buffer8 >> RUNID; throw 504;} catch (...){RUNID=1;}
	// try{char* LTANAMEchar;LTANAMEchar=((TLeafC *) config->GetBranch("LTANAME")->GetLeaf("string"))->GetValueString();istringstream buffer9(LTANAMEchar); buffer9 >> LTANAME; throw 05;} catch (...){LTANAME=1;}
}


int main(int argc, char* argv[]){

	TStopwatch t; // time counter
	t.Start();

	std::vector<string> inFileList;
	for(int i=optind; i<argc; ++i){
		inFileList.push_back(argv[i]);
	}
    std::string outfilename="outfile.root";
    TFile *outfile = new TFile(("files/"+outfilename).c_str(),"RECREATE"); //output file


    vector<vector <float>> rates(4);
    vector<vector <float>> ratesErr(4);
    vector<vector <int>> entries(4);
    TCanvas * canvas   = new TCanvas("Plots","Plots",0,0,2000,1000); //coord top lef, size x, sizey
    canvas->Print((outfilename+"_fitpeak.pdf[").c_str());
    for (size_t iFile = 0; iFile < inFileList.size(); iFile++)
	{
        // extract TTree from TFiles and read configuration
		TFile *file = TFile::Open(inFileList[iFile].c_str());
		cout << "Infile: " << inFileList[iFile] << endl;
		TTree *calPixTree = (TTree*)file->Get("calPixTree");
		readConfig(file);
        string cuts="!(mask & (4+16+128+256+512+1024+4096+8192+16384)) && min(distance,edgedist)>=60"; //"!(mask & (1+4+16+128+256+512+1024+4096)) && distance1e>=20 && min(distance,edgedist)>=" for 2e-
        // cuts = "edgedist >=400";
        cout <<cuts<<endl;
        calPixTree->SetAlias("edgedist",("min(min(x-"+std::to_string(CCDNPRES+1)+", "+std::to_string((CCDNCOL/2+CCDNPRES)/NBINCOL)+"-x), min(y-1, "+std::to_string(min(NROW,CCDNROW/(2*NBINROW)))+"-y))").c_str());


        //fit 1e peak

        std:string funcformulas="";
        for (size_t i = 0; i < 2; i++)
        {
            funcformulas=(funcformulas+"+[0]*(TMath::Gaus(x,[1]+"+std::to_string(i)+"*[2],[2]*[3],1)*TMath::Poisson("+std::to_string(i)+",[4]))").c_str();
        }
        TF1 * fitfunc = new TF1("fitfunc",funcformulas.c_str());
        std::string binning="(250,-1.3,4.0)";   
        float binwidth = (4.0-(-1.3))/250;
        
        gStyle->SetOptFit(1);
        gStyle->SetStatW(0.25); gStyle->SetStatH(0.25);
        canvas->cd();
        canvas->Divide(2,2);
        outfile->cd();
        for (size_t hdu = 1; hdu < 5; hdu++)
        {
            canvas->cd(hdu);
            gPad->SetLogy(1);
            calPixTree->Draw(("ePix>>hname"+std::to_string(hdu)+binning).c_str(),(cuts+" && ohdu=="+std::to_string(hdu)).c_str());
            TH1F * h0 = new TH1F();
            h0 = (TH1F*)outfile->Get(("hname"+std::to_string(hdu)).c_str());
            float aduZero=0; float aduPerElectron=1; float noise=0.14;
            fitfunc->SetParameters(h0->GetEntries()*binwidth,aduZero,aduPerElectron,noise,0.00001);
            h0->Fit(fitfunc,"QSL","",aduZero-1*aduPerElectron,aduZero+1.5);
            TF1 *fit = h0->GetFunction("fitfunc");
            rates[hdu-1].push_back(fit->GetParameter(4));
            // rates[hdu-1].push_back(calPixTree->GetEntries("ePix>0.5")/h0->GetEntries());
            ratesErr[hdu-1].push_back(fit->GetParError(4));
            entries[hdu-1].push_back(h0->GetEntries());

        }
        canvas->cd();
        canvas->Print((outfilename+"_fitpeak.pdf").c_str());
        canvas->Clear();
        


    }
    canvas->Print((outfilename+"_fitpeak.pdf]").c_str());
    outfile->Write();


    t.Stop();
    t.Print();

    cout << "{";
    for (size_t i = 0; i < rates[0].size(); i++) cout << rates[0][i] << ",";
    cout << "}";

    cout << "{";
    for (size_t i = 0; i < ratesErr[0].size(); i++) cout << ratesErr[0][i]<< ",";
    cout << "}";

    cout << "{";
    for (size_t i = 0; i < entries[0].size(); i++) cout << entries[0][i]<< ",";
    cout << "}";

    // rates[0]={0.000228155,0.00017877,0.000193245,0.000165712,0.000209894,0.000143151,0.000249132,0.000218208,0.000178267,0.000191552,0.000226433,0.000179477,0.000200696,0.000304744,0.000199288,0.000271985,0.000209547,0.00021034,0.000247789,0.000193612,0.000215789,0.000265461,0.000192279,0.000182461,0.000194932,0.000220191,0.000233324,0.000176973};
    // ratesErr[0]={4.13138e-05,4.05169e-05,3.99116e-05,3.05841e-05,3.73217e-05,3.48572e-05,3.98702e-05,4.04226e-05,4.11778e-05,4.01259e-05,4.78205e-05,3.50752e-05,3.85352e-05,5.47189e-05,4.46219e-05,5.65374e-05,4.77826e-05,4.47842e-05,4.46942e-05,3.20295e-05,3.8934e-05,4.83864e-05,3.81822e-05,3.90017e-05,4.93141e-05,4.1533e-05,4.74059e-05,3.49431e-05};
    // entries[0]={130486,114595,130804,182004,156562,125553,158455,133840,109309,122010,106470,155327,148088,104280,109216,95363,93920,97935,131634,181639,148597,113361,127489,117170,82287,134342,101212,169512};

    cout << rates[0].size() << endl;
    cout << ratesErr[0].size() << endl;
    cout << entries[0].size() << endl;
    

    TCanvas * canvas2   = new TCanvas("Plots","Plots",0,0,2000,1000); //coord top lef, size x, sizey
    canvas2->Print((outfilename+"_test.pdf[").c_str());
    canvas2->cd();
    // canvas2->Divide(5,6);
    canvas2->Divide(2,2);
    gRandom->SetSeed(0);
    outfile->cd();
    float N=CCDNROW*CCDNCOL/4;
    cout << "total pixels" << N << endl;


    // KS test
    vector <double> pvalues;
    // create "theoretical" distribution
    int nsim =1000;
    double ladder[nsim][3];
    for (size_t i = 0; i < nsim; i++)
    {
        ladder[i][0]=i/float(nsim);
        ladder[i][1]=i/float(nsim);
        ladder[i][2]=1; // 0 if data 1 if theory
    }


    /////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////
    for (size_t hdu = 1; hdu < 5; hdu++)
    {
        pvalues.clear();
        TH1F * h1 = new TH1F("h1","h1",10,0,1);
        
        canvas2->cd(hdu);
        float U = (N)*weighted_average(ratesErr[hdu-1],rates[hdu-1],rates[hdu-1].size());
        cout << "U: " << U << endl;
        // Int_t NSIMS=10000000;
        
        for (size_t j = 0; j < rates[hdu-1].size(); j++)
        {
            // float nMean = accumulate(entries[hdu-1].begin(), entries[hdu-1].end(), 0.0)/entries[hdu-1].size();

            // float n=0;
            // // cout << "nMean: " << entries[hdu-1][j] << endl;
            // for (size_t k= 0; k < NSIMS; k++)
            // {
            //     // float A = gRandom->Poisson(U);
            //     // float a = gRandom->Gaus(U,pow(U,0.5)*(N/entries[hdu-1][j]))/(N);
            //     // float A = gRandom->Gaus(U,pow(U,0.5)*(N/entries[hdu-1][j]));
            //     float A = gRandom->Gaus(U,pow(U*(N/entries[hdu-1][j]),0.5));
                
            //     if (A>rates[hdu-1][j]*(N))
            //     {
            //         n++;
            //     }

            // }
            // cout << n/NSIMS << endl;

            TF1 fGaus_cdf("f","ROOT::Math::normal_cdf_c(x,[0],[1])",0,1000); // [0,vObs[i]*2+10] is the range
            // fGaus_cdf.SetParameter(0,pow(U*(N/entries[hdu-1][j]),0.5));
            // fGaus_cdf.SetParameter(1,U);
            // double p = fGaus_cdf(rates[hdu-1][j]*N);

            fGaus_cdf.SetParameter(0,pow(entries[hdu-1][j]*U/N,0.5));
            fGaus_cdf.SetParameter(1,entries[hdu-1][j]*U/N);
            double p = fGaus_cdf(rates[hdu-1][j]*entries[hdu-1][j]);


            cout << "File: " << inFileList[j] << endl;
            cout << p << endl;
            cout << "------------" << endl;
            // if (p<0.002 || p>1-0.002)
            // {
            //     continue;
            // }
            // h1->Fill(n/NSIMS);
            // pvalues.push_back(n/NSIMS);
            h1->Fill(p);
            pvalues.push_back(p);
            

            
            
        }

        //// create experimental ladder
        sort(pvalues.begin(), pvalues.end());
        double ladderexp[int(pvalues.size())][3];
        for (size_t j = 0; j < pvalues.size(); j++)
        {
            ladderexp[j][0]=pvalues[j];
            ladderexp[j][1]=j/float(pvalues.size());
            ladderexp[j][2]=0; // 0 if data 1 if theory
        }
        h1->Draw();
        

        ////////////////////////
        // float simStdDev=h1->GetStdDev();
        // vector <double> measuredStdDevVector;
        // for (size_t j = 0; j < rates[hdu-1].size(); j++)
        // {
        //     TLine *l = new TLine(rates[hdu-1][j]*(3072*512),100,rates[hdu-1][j]*(3072*512),300);
        //     measuredStdDevVector.push_back(rates[hdu-1][j]*(3072*512));
        //     l->Draw();
        // }
        
        // float measuredStdDev = stddev(measuredStdDevVector);


        
        double escalera[nsim+pvalues.size()][3];
        concatenate_doubles(escalera, ladder,ladderexp,nsim,pvalues.size());
        qsort(escalera, ARRSIZE(escalera), sizeof(*escalera), compare);

        double largest_diff=0;
        compute_largest_diff(escalera, pvalues, largest_diff);

        // cout << largest_diff << endl;
        // cout << "--------------------------------" << endl;

        double real_largest_diff=largest_diff;
        double counter=0;
        auto rndm = TRandom3(0);
	    rndm.SetSeed(0);
        for (size_t i = 0; i < 10000; i++)
        {
            vector <double> datasim;
            for (size_t j = 0; j < pvalues.size(); j++)
            {
                datasim.push_back(rndm.Uniform(0,1));
            }
            sort(datasim.begin(), datasim.end());
            double laddersim[int(pvalues.size())][3];
            for (size_t j = 0; j < datasim.size(); j++)
            {
                laddersim[j][0]=datasim[j];
                laddersim[j][1]=j/float(datasim.size());
                laddersim[j][2]=0; // 0 if data 1 if theory
            }
            double escalera[nsim+pvalues.size()][3];
            concatenate_doubles(escalera, ladder,laddersim,nsim,pvalues.size());
            qsort(escalera, ARRSIZE(escalera), sizeof(*escalera), compare);

            double largest_diff=0;
            compute_largest_diff(escalera,pvalues, largest_diff);

            // cout << largest_diff << endl;
            // cout << "--------------------------------" << endl;  
            if (largest_diff>real_largest_diff)
            {
                counter++;
            }
              
        }

        char tv_KS[40];
        sprintf(tv_KS,  "p-value = %3.2f",counter/10000);
        TPaveLabel *px = new TPaveLabel(0.62,0.15,0.83,0.24,tv_KS,"brNDC");
        px->SetBorderSize(0);px->SetTextFont(42); px->SetTextSize(0.5);px->SetTextAlign(22);px->SetTextColor(1); 
        px->Draw();
        cout << counter/100 << endl;
        
        



    }
    
    canvas2->cd();
    canvas2->Print((outfilename+"_test.pdf").c_str());
    canvas2->Clear();
    canvas2->Print((outfilename+"_test.pdf]").c_str());
    outfile->Write();
    
    
    



    /////////////////////////////////////// Likelihood test ////////////////////////////////////////


    // //initialize
    // gRandom->SetSeed(0);
    // for (size_t hdu = 1; hdu < 5; hdu++)
    // {
    //     vector <float> mutemp;
    //     for (size_t i = 0; i < rates[hdu-1].size(); i++) mutemp.push_back(rates[hdu-1][i]*entries[hdu-1][i]/accumulate(entries[hdu-1].begin(), entries[hdu-1].end(), 0.0));
        
    //     float mu = accumulate( mutemp.begin(), mutemp.end(), 0.0);  
    //     // float mu = accumulate( rates[hdu-1].begin(), rates[hdu-1].end(), 0.0)/rates[hdu-1].size();  
    //     Int_t NSIMS=10000;
    //     vector<vector<Int_t>> mat;
    //     mat.resize(NSIMS, std::vector<Int_t>(rates[hdu-1].size()));

    //     // simulate events
    //     for (size_t i = 0; i < rates[hdu-1].size(); i++)
    //     {
    //         for (size_t j = 0; j < NSIMS; j++)
    //         {
    //             int a = gRandom->Poisson(mu*entries[hdu-1][i]);
    //             // cout << a << endl;
    //             mat[j][i]=a; // each column is the same bin
    //         }
    //         // cout << "************************" << endl;
    //     }

    //     // calculate likelihoods
    //     vector<Float_t> Lj;
    //     vector<Float_t> pi;
    //     for (vector<Int_t> row: mat) {
    //         pi.clear();
    //         int i=0;
    //         for (Int_t value : row){
    //             // DO NOT DELETE THIS COMMENT BELOW!!
    //             // pi.push_back(-log(exp(-mu*entries[hdu-1][i])*pow(mu*entries[hdu-1][i],value)/tgamma(value+1))); 
    //             pi.push_back(mu*entries[hdu-1][i]-value*log(mu*entries[hdu-1][i])+lgamma(value+1));
    //             // cout << "pi: " <<  mu*entries[hdu-1][i]-value*log(mu*entries[hdu-1][i])+lgamma(value+1) << endl;
    //             // cout << "value: "<< value <<endl;
    //             i++;
    //         }
    //         // cout << " | "<<endl;
    //         Lj.push_back(std::accumulate(pi.begin(), pi.end(), 0.00000));
    //         // cout << std::accumulate(pi.begin(), pi.end(), 0.00000) << endl;
    //         // cout << "--------------------------------" << endl;
    //     }
    //     cout << "Mean value: "<< std::accumulate(Lj.begin(), Lj.end(), 0.00000)/Lj.size() << endl;
    //     double max = *max_element(Lj.begin(), Lj.end());
    //     cout<<"Max value : "<<max<<endl;

    //     // calculate measured likelihood
        
    //     pi.clear();
    //     for (size_t i = 0; i < entries[hdu-1].size(); i++)
    //     {

    //     // pi.push_back(-log(exp(-rate[i]*entries[hdu-1][i])*pow(rate[i]*entries[hdu-1][i],rates[hdu-1][i]*entries[hdu-1][i])/tgamma(rates[hdu-1][i]*entries[hdu-1][i]+1)));
    //     // pi.push_back(-log(exp(-mu*entries[hdu-1][i])*pow(mu*entries[hdu-1][i],rates[hdu-1][i]*entries[hdu-1][i])/tgamma(rates[hdu-1][i]*entries[hdu-1][i]+1)));
    //     pi.push_back(mu*entries[hdu-1][i]-rates[hdu-1][i]*entries[hdu-1][i]*log(mu*entries[hdu-1][i])+lgamma(rates[hdu-1][i]*entries[hdu-1][i]+1));
    //     // cout << mu*entries[hdu-1][i]-rates[hdu-1][i]*entries[hdu-1][i]*log(mu*entries[hdu-1][i])+lgamma(rates[hdu-1][i]*entries[hdu-1][i]+1) << endl;
    //     // cout << mu*entries[hdu-1][i] << endl;
    //     // cout << rates[hdu-1][i]*entries[hdu-1][i] << endl;
    //     // cout << -rates[hdu-1][i]*entries[hdu-1][i]*log(mu*entries[hdu-1][i]) << endl;
    //     // cout << lgamma(rates[hdu-1][i]*entries[hdu-1][i]+1) << endl;
    //     // cout << "%%%%%%" << endl;

    //     // cout << rate[i]*entries[hdu-1][i] << " | " <<" | " << pi[i] << endl<< endl;
    //     // cout << pi[i] << " | " << rate[i]*entries[hdu-1][i] << " | "<< rates[hdu-1][i]*entries[hdu-1][i] << " | " << i <<  endl<< endl;
    //     }


    //     Double_t pi_value= std::accumulate(pi.begin(), pi.end(), 0.00000);
    //     Float_t nLikelihood=0;
    //     cout << "my p-value: " <<  std::accumulate(pi.begin(), pi.end(), 0.00000) << endl;
    //     for(float i : Lj) if(i>=pi_value){nLikelihood++;}
    //     cout << "p-value Likelihood = " << nLikelihood/NSIMS << endl;
    //     cout << nLikelihood << " | "  << NSIMS << endl; 
    // }
    
    return 0;
}