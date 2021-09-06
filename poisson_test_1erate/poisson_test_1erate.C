#include "TPaveLabel.h"
#include "TLegend.h"
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

// concatenate double variables
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

// calculate largest distance between to cumulative distribution to perform K-S test
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

// weighted_average of a vector //
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

// Standard deviation of a vector//
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


void count2e(std::vector<string> &inFileList,size_t &hdu, vector<vector<float>> &two_e_1pix, vector<vector<float>> &two_e_h, vector<vector<float>> &two_e_v, vector<vector<float>> &two_e_d, vector<vector<float>> & three_e)
{
    cout << "Counting 2e- events" << endl;
    for (size_t iFile = 0; iFile < inFileList.size(); iFile++)
	{
        // extract TTree from TFiles and read configuration
		TFile *file = TFile::Open(inFileList[iFile].c_str());
		cout << "Infile: " << inFileList[iFile] << endl;
		TTree *hitSumm = (TTree*)file->Get("hitSumm");
        // string cuts2e="!(flag & (4+16+64+128+256+512+1024+4096+8192+16384)) && distance1e>=30 && distance>=30";
        string cuts2e="!(flag & (4+16+64+128+256+512+1024+4096+8192+16384)) && distance>=30";

        two_e_1pix[hdu-1].push_back(hitSumm->GetEntries((cuts2e+" && e==2 && n==1 && ohdu=="+std::to_string(hdu)).c_str()));
        two_e_h[hdu-1].push_back(hitSumm->GetEntries((cuts2e+" && e==2 && n==2 && yVar==0 && ohdu=="+std::to_string(hdu)).c_str()));
        two_e_v[hdu-1].push_back(hitSumm->GetEntries((cuts2e+" && e==2 && n==2 && xVar==0 && ohdu=="+std::to_string(hdu)).c_str()));
        two_e_d[hdu-1].push_back(hitSumm->GetEntries((cuts2e+" && e==2 && n==2 && xVar!=0 && yVar!=0 && ohdu=="+std::to_string(hdu)).c_str()));
        three_e[hdu-1].push_back(hitSumm->GetEntries((cuts2e+" && e==3 && ohdu=="+std::to_string(hdu)).c_str()));

    }
    
    cout << two_e_1pix[hdu-1].size()<< endl;
    cout << two_e_h[hdu-1].size()<< endl;
    cout << two_e_v[hdu-1].size()<< endl;
    cout << two_e_d[hdu-1].size()<< endl;
    cout << three_e[hdu-1].size()<< endl;
    return;

}

int main(int argc, char* argv[]){

	TStopwatch t; // time counter
	t.Start();

	std::vector<string> inFileList;
    std::vector<string> inFileList2;
	for(int i=optind; i<argc-3; ++i){
		inFileList.push_back(argv[i]);
        std::size_t pos = inFileList.back().find("run"); // position of "run" in argv[i]
        inFileList2.push_back(inFileList.back().substr(pos));   
        cout << argv[i] << endl;
	}
    std::string outfilename= argv[argc-3];
    cout << outfilename << endl;
    TFile *outfile = new TFile(("files/"+outfilename+".root").c_str(),"RECREATE"); //output file

    float pvalue_thr = atof(argv[argc-2]);


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
        string cuts="!(mask & (1+4+16+64+128+256+512+1024+4096+8192+16384)) && distance>=60"; //"!(mask & (1+4+16+128+256+512+1024+4096)) && distance1e>=20 && min(distance,edgedist)>=" for 2e-
        // cuts = "edgedist >=450";
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

            cout << "{";
            for (size_t i = 0; i < rates[hdu-1].size(); i++) cout << rates[hdu-1][i] << ",";
            cout << "}";

            cout << "{";
            for (size_t i = 0; i < ratesErr[hdu-1].size(); i++) cout << ratesErr[hdu-1][i]<< ",";
            cout << "}";

            cout << "{";
            for (size_t i = 0; i < entries[hdu-1].size(); i++) cout << entries[hdu-1][i]<< ",";
            cout << "}";

        }
        canvas->cd();
        canvas->Print((outfilename+"_fitpeak.pdf").c_str());
        canvas->Clear();
        


    }
    canvas->Print((outfilename+"_fitpeak.pdf]").c_str());
    outfile->Write();


    t.Stop();
    t.Print();

    
    // std::string ratesstring = argv[argc-1];
    // if (ratesstring=="0module1")
    // {
    //     rates[0]={9.66523e-05,0.000102099,0.000105092,7.63878e-05,9.61424e-05,8.84042e-05,9.3425e-05,8.8555e-05,8.77373e-05,7.88392e-05,8.88468e-05,0.00010336,0.000103503,7.79742e-05,7.99184e-05,9.3088e-05,9.50554e-05,6.65978e-05,6.83975e-05,8.38009e-05,6.51527e-05};ratesErr[0]={1.19056e-05,1.27353e-05,1.18324e-05,9.54987e-06,1.12907e-05,1.05249e-05,1.08726e-05,1.05636e-05,1.10762e-05,9.82206e-06,1.00974e-05,1.15456e-05,1.1584e-05,9.95879e-06,1.06782e-05,1.2154e-05,1.16523e-05,9.89047e-06,9.8002e-06,1.08476e-05,9.20504e-06};entries[0]={709382,724170,798518,868788,793525,832600,759775,853214,748879,846669,864873,820732,830938,854372,746518,705066,793125,805010,745526,785626,827739};rates[1]={6.85062e-05,6.21974e-05,0.000109854,0.000111951,8.74784e-05,8.8527e-05,9.97559e-05,8.76845e-05,5.12124e-05,7.24843e-05,6.72173e-05,6.04912e-05,8.33669e-05,9.53638e-05,8.74316e-05,9.49465e-05,8.90224e-05,5.36285e-05,8.5495e-05,8.58577e-05,6.90976e-05};ratesErr[1]={1.05627e-05,9.23383e-06,1.2947e-05,1.22836e-05,1.12647e-05,1.04786e-05,1.17979e-05,1.16404e-05,7.79632e-06,9.93806e-06,9.89495e-06,9.1292e-06,1.12346e-05,1.2098e-05,1.00441e-05,1.20312e-05,1.15525e-05,8.07647e-06,1.15113e-05,1.11984e-05,9.35273e-06};entries[1]={682472,779871,702372,711817,737174,870415,743667,705240,884848,746095,756782,793694,731577,670658,915423,679107,725201,884853,662383,736922,816243};rates[2]={8.40747e-05,7.39293e-05,8.55842e-05,9.34596e-05,8.39588e-05,9.49696e-05,8.96245e-05,6.87097e-05,9.48092e-05,9.00961e-05,8.33049e-05,8.18049e-05,7.6529e-05,7.2515e-05,8.03264e-05,7.63327e-05,7.83967e-05,9.50188e-05,6.23664e-05,6.16052e-05,6.16225e-05};ratesErr[2]={9.76396e-06,9.00281e-06,9.45314e-06,9.92e-06,9.73434e-06,1.00825e-05,9.8926e-06,8.7123e-06,9.94841e-06,9.98423e-06,1.00002e-05,9.94857e-06,9.55362e-06,9.32421e-06,1.00945e-05,1.01007e-05,9.58444e-06,1.01973e-05,8.49369e-06,8.31467e-06,8.38513e-06};entries[2]={909666,937679,954069,965049,907296,948986,940823,949403,968296,922681,834820,861106,860845,876073,807550,761475,897396,897610,896637,906852,905434};rates[3]={8.17527e-05,9.77451e-05,9.91465e-05,9.3067e-05,9.97123e-05,9.83065e-05,0.000107162,0.000103253,6.71519e-05,9.67408e-05,9.50447e-05,8.98694e-05,9.04098e-05,6.77779e-05,8.49418e-05,8.77207e-05,7.9233e-05,8.24396e-05,8.72502e-05,0.000115985,7.46644e-05};ratesErr[3]={9.48092e-06,1.14553e-05,1.14189e-05,1.08573e-05,1.05019e-05,1.11729e-05,1.16748e-05,1.1039e-05,9.40267e-06,1.06952e-05,1.09844e-05,1.00271e-05,1.04319e-05,9.77349e-06,1.01863e-05,1.03863e-05,9.46004e-06,1.09554e-05,1.06259e-05,1.29691e-05,9.88866e-06};entries[3]={939576,885986,838195,835944,956660,870869,855098,914479,800442,833408,895982,961099,897637,777770,875028,875147,915769,718995,828781,780081,800876};
    //     cout << rates[0].size() << rates[1].size() << rates[2].size() << rates[3].size() << endl;
    //     cout << ratesErr[0].size() << ratesErr[1].size() << ratesErr[2].size() << ratesErr[3].size() << endl;
    //     cout << entries[0].size() << entries[1].size() << entries[2].size() << entries[3].size() << endl;
    // }
    // if (ratesstring=="0module2")
    // {
    //     rates[0]={0.000731692,0.000304154,0.000617719,0.000643935,0.000797164,0.000622017,0.000861024,0.00291055,0.00266388,0.00146138,0.000907955,0.00100952,0.000831136,0.00132528,0.00109517,0.000934319,0.00129967,0.00136757,0.00192638,0.000967035,0.00117304};ratesErr[0]={3.89762e-05,2.51446e-05,3.55228e-05,3.54736e-05,4.249e-05,3.3237e-05,4.42384e-05,7.99634e-05,7.5971e-05,5.10695e-05,3.98543e-05,4.55214e-05,4.504e-05,5.50673e-05,4.59183e-05,4.50345e-05,5.06348e-05,5.23325e-05,6.32536e-05,4.45459e-05,5.05732e-05};entries[0]={501404,489614,510674,536883,461739,586787,509306,471379,474130,546370,599927,508489,428132,448091,550038,486620,518546,493360,497618,518607,477370};rates[1]={0.00251763,0.00106838,0.00214487,0.00201567,0.00236125,0.00201266,0.00254315,0.00253518,0.00276346,0.00232041,0.00206204,0.00196202,0.00177333,0.00226831,0.00229296,0.00223439,0.00263305,0.0029983,0.00253737,0.00232288,0.00363428};ratesErr[1]={7.55325e-05,5.04256e-05,6.93937e-05,6.24205e-05,6.79971e-05,6.36385e-05,8.44247e-05,0.000120233,0.000122091,9.95905e-05,7.87539e-05,7.60667e-05,7.39685e-05,9.30695e-05,8.80763e-05,0.000109403,0.00010156,0.000119725,0.000107268,8.95692e-05,0.000134607};entries[1]={458998,453635,461000,536927,533107,515134,367829,181650,190953,244282,343756,352362,335923,272640,307247,195763,264552,215803,226836,299008,207040};rates[2]={0.000338463,0.000168919,0.000284978,0.000276768,0.000342885,0.000297119,0.00037972,0.000423709,0.000566098,0.000382309,0.00031467,0.00034013,0.00028213,0.000444519,0.000461426,0.000455271,0.00056395,0.000513521,0.00049021,0.000397011,0.000516515};ratesErr[2]={2.00004e-05,1.49574e-05,1.86274e-05,1.98261e-05,2.06538e-05,1.99413e-05,2.07727e-05,2.29631e-05,2.63845e-05,2.17373e-05,1.88869e-05,1.95203e-05,1.90016e-05,2.35284e-05,2.34549e-05,2.2455e-05,2.56438e-05,2.47213e-05,2.52931e-05,2.28636e-05,2.57632e-05};entries[2]={920700,774999,904659,770377,896811,852823,953397,863020,824430,883478,946855,971134,887655,860483,866328,954550,939548,865431,823030,841330,825559};rates[3]={0.000180426,0.000107069,0.000174988,0.000160201,0.000164987,0.000175134,0.000214981,0.000474433,0.0004787,0.0002921,0.000207974,0.000231168,0.000173683,0.000202146,0.000221082,0.000218505,0.000243547,0.000274547,0.000370134,0.000212583,0.00025054};ratesErr[3]={1.53369e-05,1.18259e-05,1.36439e-05,1.4186e-05,1.46852e-05,1.39353e-05,1.56081e-05,2.46278e-05,2.28862e-05,1.99742e-05,1.60254e-05,1.71167e-05,1.41681e-05,1.55972e-05,1.66897e-05,1.66202e-05,1.74174e-05,1.89173e-05,2.15744e-05,1.62634e-05,1.78313e-05};entries[3]={833501,822308,981200,890975,839738,924311,922873,854920,960902,761504,840778,826928,943345,893948,853265,803570,862889,797022,830997,843253,857819};
    //     cout << rates[0].size() << rates[1].size() << rates[2].size() << rates[3].size() << endl;
    //     cout << ratesErr[0].size() << ratesErr[1].size() << ratesErr[2].size() << ratesErr[3].size() << endl;
    //     cout << entries[0].size() << entries[1].size() << entries[2].size() << entries[3].size() << endl;
    // }
    // if (ratesstring=="72000module1")
    // {
    //     rates[0]={0.000228579,0.000179047,0.000193622,0.000165954,0.000197702,0.000143351,0.000249652,0.00021115,0.000178755,0.000191829,0.000226818,0.000179752,0.00020113,0.000305445,0.000199592,0.000272559,0.00020986,0.000210728,0.000248337,0.000193995,0.000216189,0.000249954,0.000165594,0.000182695,0.000195272,0.000220601,0.000233921,0.000177287};ratesErr[0]={4.13866e-05,4.06211e-05,4.00809e-05,3.14537e-05,3.6504e-05,3.49656e-05,3.99486e-05,3.95751e-05,4.12911e-05,4.04768e-05,4.74722e-05,3.51718e-05,3.85789e-05,5.45095e-05,4.31532e-05,5.51712e-05,4.80024e-05,4.49054e-05,4.46849e-05,3.27182e-05,3.90306e-05,4.75972e-05,3.67153e-05,4.11356e-05,4.93105e-05,4.12888e-05,4.88206e-05,3.34353e-05};entries[0]={130256,114418,130576,181758,156284,125386,158135,133610,109030,121824,106295,155109,147792,104033,109055,95167,93770,97759,131370,181350,148344,113099,127289,117003,82152,134100,100970,169233};rates[1]={0.000209933,0.000116424,0.000191586,0.000192258,0.00014837,0.000201112,0.000234192,0.000208731,0.000230657,0.000166679,0.000224145,0.000159902,0.000148622,0.000198193,0.000264587,0.000236303,0.00018548,0.000137831,0.0001827,0.000237549,0.000106162,0.000178709,0.000228108,0.00036547,0.000176002,0.000280234,0.000223205,0.000115485};ratesErr[1]={4.75721e-05,5.26799e-05,4.37211e-05,4.89756e-05,4.15301e-05,3.66177e-05,6.52713e-05,4.58735e-05,4.59547e-05,4.18171e-05,5.46524e-05,4.67521e-05,4.952e-05,4.50693e-05,5.18957e-05,4.67921e-05,5.29546e-05,4.37867e-05,7.4274e-05,6.1105e-05,3.31512e-05,5.35842e-05,6.50674e-05,6.47474e-05,4.68924e-05,6.3043e-05,5.1047e-05,4.8061e-05};entries[1]={94205,46187,115876,78379,91562,155904,57582,104899,115513,96502,76440,82107,61014,100526,102606,106652,73211,84738,33478,60159,97629,69537,56942,91315,80915,74203,94987,54518};rates[2]={0.000126219,0.000170343,0.000184397,0.000156949,0.000145633,0.000182651,0.00018236,0.000216848,0.000219266,0.000150108,0.000179172,0.000191016,0.000159018,0.000177806,0.000174358,0.000323248,0.000147057,0.000157568,0.000130523,0.000124382,0.000175841,0.00016326,0.000130363,0.000204274,0.000141468,0.000146735,0.000222797,0.000161304};ratesErr[2]={2.15582e-05,2.42397e-05,2.88616e-05,2.19298e-05,2.38508e-05,2.60674e-05,3.22952e-05,3.52182e-05,3.02077e-05,2.44527e-05,2.71503e-05,3.14132e-05,2.85216e-05,3.03192e-05,2.62854e-05,5.30793e-05,2.31316e-05,2.42817e-05,2.56083e-05,2.30297e-05,3.39907e-05,2.97789e-05,2.46056e-05,3.27797e-05,3.13824e-05,2.97633e-05,3.21897e-05,2.86076e-05};entries[2]={272294,295236,228436,325623,258204,270093,175549,176220,246326,252856,245096,194320,199301,196631,259098,121982,281540,270425,206324,244093,153666,186206,211473,192257,143517,163408,215925,198669};rates[3]={0.00016779,0.000187638,0.000168466,0.000145383,0.000159144,0.000213328,0.000188819,0.000161626,0.000216737,0.000176359,0.000213376,0.000190369,0.000184284,0.000159161,0.000109114,0.000152272,0.000190585,0.000176879,0.000190505,0.000217307,0.000146523,0.000149203,0.000144642,0.000202696,0.000164549,0.000230057,0.000259389,0.000215668};ratesErr[3]={4.01534e-05,3.06355e-05,2.47209e-05,2.39879e-05,2.37227e-05,3.35677e-05,2.96645e-05,2.67256e-05,3.8152e-05,2.33978e-05,4.05423e-05,3.0433e-05,2.91657e-05,2.6664e-05,2.22417e-05,3.00571e-05,2.71612e-05,2.5858e-05,3.20009e-05,3.09331e-05,3.15565e-05,3.36862e-05,2.98439e-05,4.00828e-05,3.26722e-05,3.47041e-05,3.78052e-05,3.65851e-05};entries[3]={108678,210962,297330,282928,312096,201435,241304,267660,162732,347641,135996,233859,227316,247608,229657,178365,273881,283075,193882,244117,157764,134039,169047,132901,176303,218896,195414,174632};
    //     cout << rates[0].size() << rates[1].size() << rates[2].size() << rates[3].size() << endl;
    //     cout << ratesErr[0].size() << ratesErr[1].size() << ratesErr[2].size() << ratesErr[3].size() << endl;
    //     cout << entries[0].size() << entries[1].size() << entries[2].size() << entries[3].size() << endl;
    // }
    // if (ratesstring=="72000module2")
    // {
    //     rates[0]={0.0013553,0.00106247,0.00103645,0.000848389,0.00145374,0.00115667,0.00122894,0.00130336,0.000982585,0.000870627,0.00121538,0.000818379,0.00121264,0.00168207,0.00120917,0.00168063,0.00130725,0.00105595,0.00208799,0.00259832,0.00496777,0.00183778,0.00271529,0.00171347,0.00138352,0.00183588,0.0012654,0.00117678};ratesErr[0]={0.000117909,0.000106618,0.00015332,8.89706e-05,0.000155558,9.55556e-05,0.000123876,0.000117643,9.94517e-05,8.92806e-05,0.000141463,0.000100988,0.000137147,0.000133749,9.61891e-05,0.000132247,9.53045e-05,9.53464e-05,0.000185413,0.000155083,0.000244386,0.000138844,0.00015225,0.000135812,9.49006e-05,0.000146229,0.000126808,0.000121233};entries[0]={93311,94739,44412,112839,55910,129096,83163,94593,109141,114436,57242,84319,64422,93862,138085,86623,137555,123857,62686,110955,88193,96662,121562,81214,148130,80984,75980,87534};rates[1]={0.00263365,0.00290457,0.00232432,0.00223163,0.00305539,0.00250489,0.00310662,0.00286135,0.00189457,0.00197596,0.00332114,0.00252464,0.00394022,0.0021506,0.00248203,0.00279537,0.00264481,0.00261744,0.00275821,0.00346628,0.00298766,0.00247685,0.00264595,0.0026089,0.00182169,0.00273035,0.00194867,0.00210098};ratesErr[1]={0.000187649,0.000192131,0.000146504,0.000142422,0.000174906,0.00017042,0.000173928,0.000250922,0.000183448,0.000147437,0.000234076,0.000205367,0.000269811,0.000166729,0.000194582,0.000297916,0.000213993,0.000219616,0.000218993,0.000347718,0.000435697,0.000215843,0.000256268,0.000275715,0.000166723,0.000312827,0.000207321,0.000233367};entries[1]={76515,81663,113720,114130,103039,88723,105846,46690,58291,95491,62468,62244,54640,60915,67367,28865,59260,55550,59959,29670,15859,55511,42106,35728,67983,29157,46749,39497};rates[2]={0.000415344,0.000330272,0.000386517,0.000442318,0.00044365,0.000480692,0.000368483,0.000331892,0.000418223,0.00037482,0.000365955,0.000254811,0.000504142,0.000548838,0.000524486,0.000515904,0.000483696,0.000341284,0.000615797,0.000690408,0.000826152,0.000580063,0.00060321,0.00050433,0.00048416,0.000411541,0.00046217,0.000400315};ratesErr[2]={5.2094e-05,4.01268e-05,4.39856e-05,4.94031e-05,4.81842e-05,4.56025e-05,5.37471e-05,4.19196e-05,4.96134e-05,4.68457e-05,4.1441e-05,3.84453e-05,5.76378e-05,5.66415e-05,4.84529e-05,4.96092e-05,4.68909e-05,3.66828e-05,5.96949e-05,5.97569e-05,6.09296e-05,5.6823e-05,5.79623e-05,6.09485e-05,5.48535e-05,4.87267e-05,4.35431e-05,5.49075e-05};entries[2]={159818,221485,213359,192875,217822,230854,152007,213431,180389,189115,239873,179505,164280,176014,238937,227342,224678,278496,190978,213287,247794,193612,182517,137877,183549,188902,253403,148283};rates[3]={0.000203162,0.00023745,0.000245413,0.000240931,0.00025184,0.000312805,0.000322574,0.000460223,0.00022776,0.00016841,0.000276532,0.000222184,0.000304988,0.000437469,0.000287708,0.000346984,0.000365469,0.000300328,0.00033954,0.000799324,0.00118132,0.000413025,0.000848259,0.000449428,0.000347721,0.000355954,0.000318347,0.00026245};ratesErr[3]={3.3607e-05,2.92956e-05,3.69218e-05,3.4074e-05,3.25733e-05,3.63849e-05,4.19077e-05,4.71885e-05,3.358e-05,3.11337e-05,4.18852e-05,3.34799e-05,4.02434e-05,4.79197e-05,3.71243e-05,4.00491e-05,4.04534e-05,3.40366e-05,4.42904e-05,6.26386e-05,7.84339e-05,4.10844e-05,6.82006e-05,5.15997e-05,4.83532e-05,3.71048e-05,4.54861e-05,3.88682e-05};entries[3]={186444,294520,194064,208757,255001,254158,192539,214519,209919,187296,163573,202642,202968,206953,197973,211347,219617,273836,180612,229584,200581,253638,197814,176324,164158,279518,160185,187629};
    //     cout << rates[0].size() << rates[1].size() << rates[2].size() << rates[3].size() << endl;
    //     cout << ratesErr[0].size() << ratesErr[1].size() << ratesErr[2].size() << ratesErr[3].size() << endl;
    //     cout << entries[0].size() << entries[1].size() << entries[2].size() << entries[3].size() << endl;
    // }

    

    TCanvas * canvas2   = new TCanvas("Plots","Plots",0,0,2000,1000); //coord top lef, size x, sizey
    canvas2->Print((outfilename+"_test.pdf[").c_str());
    canvas2->cd();
    canvas2->Divide(2,2);
    gRandom->SetSeed(0);
    outfile->cd();
    float N=CCDNROW*CCDNCOL/4;
    cout << "Total pixels: " << N << endl;


    //2e counting
    vector<vector <float>> two_e_1pix(4);
    vector<vector <float>> two_e_h(4);
    vector<vector <float>> two_e_v(4);
    vector<vector <float>> two_e_d(4);
    vector<vector <float>> three_e(4);
    // vector<vector <float>> ratesErr(4);

    // KS test
    vector <double> pvalues;
    vector <double> pvaluestemp;
    std::vector<string> inFileListAux;
    // write into csv
    std::ofstream myfile;
    myfile.open ((outfilename+".csv").c_str());
    myfile << "filename,1e- rate,1e- rateErr,p-value,2e- count(1pix),2e- count(h),2e- count(v),2e- count(d),2e- count(v+d),3e- (all),\n";
    // create "theoretical" distribution
    int nsim =10000;
    double ladder[nsim][3];
    for (size_t i = 0; i < nsim; i++)
    {
        ladder[i][0]=i/float(nsim);
        ladder[i][1]=i/float(nsim);
        ladder[i][2]=1; // 0 if data 1 if theory
    }


    /////////////////////////////////////////////////////////////////////////////////
    for (size_t hdu = 1; hdu < 5; hdu++)
    {

        count2e(inFileList,hdu, two_e_1pix, two_e_h, two_e_v, two_e_d, three_e);

        inFileListAux=inFileList2;
        bool go = true;
        while (go)
        {       
            pvalues.clear(); pvaluestemp.clear();
            TH1F * h1 = new TH1F("h1","h1",10,0,1);
            
            canvas2->cd(hdu);
            float u = weighted_average(ratesErr[hdu-1],rates[hdu-1],rates[hdu-1].size());
            cout << "u: " << u << endl;
            cout << "len(u): " << rates[hdu-1].size() << endl;
            cout << "{";
            for (size_t i = 0; i < rates[hdu-1].size(); i++) cout << rates[hdu-1][i] << ",";
            cout << "}";
            
            for (size_t j = 0; j < rates[hdu-1].size(); j++)
            {

                TF1 fGaus_cdf("f","ROOT::Math::normal_cdf_c(x,[0],[1])",0,1000); // [0,vObs[i]*2+10] is the range
                gStyle->SetOptStat(1);

                fGaus_cdf.SetParameter(0,pow(entries[hdu-1][j]*u,0.5));
                fGaus_cdf.SetParameter(1,entries[hdu-1][j]*u);
                double p = fGaus_cdf(entries[hdu-1][j]*rates[hdu-1][j]);

                if (p>0.5)
                {
                    pvaluestemp.push_back(1-p);
                }
                else
                {
                    pvaluestemp.push_back(p);
                }

                h1->Fill(p);
                pvalues.push_back(p);
         
                
            }

            
            // // ROOT::Math::binomial_cdf_c (unsigned int k, double p, unsigned int n)
            // TF1 fBinomial_cdf_c("f","ROOT::Math::binomial_cdf_c(x,[0],[1])",0,20);
            // fBinomial_cdf_c.SetParameter(0,pvalue_thr);
            // fBinomial_cdf_c.SetParameter(1,pvalues.size());

            // int contador=0;
            // for (size_t j = 0; j < pvalues.size(); j++)
            // {
            //     if (pvalues[j]<pvalue_thr)
            //     {
            //         contador++;
            //     }
            // }

            // // cout << contador << endl;
            // // cout << pvalues.size()*pvalue_thr << endl;
            // // cout << 2*pow(pvalues.size()*pvalue_thr*pvalue_thr*(1-pvalue_thr),0.5) << endl;
            // // if (contador>pvalues.size()*pvalue_thr+2*pow(pvalues.size()*pvalue_thr*pvalue_thr*(1-pvalue_thr),0.5))


            // cout << contador << endl;
            // cout << fBinomial_cdf_c(contador) << endl;
            // cout << pvalues.size()*pvalue_thr << endl;
            // if (contador>pvalues.size()*pvalue_thr)
            // // if (fBinomial_cdf_c(contador)<0.0001)
            // {
            //     rates[hdu-1].erase(rates[hdu-1].begin()+std::distance(pvalues.begin(),std::min_element(pvalues.begin(), pvalues.end())));
            //     ratesErr[hdu-1].erase(ratesErr[hdu-1].begin()+std::distance(pvalues.begin(),std::min_element(pvalues.begin(), pvalues.end())));
            //     entries[hdu-1].erase(entries[hdu-1].begin()+std::distance(pvalues.begin(),std::min_element(pvalues.begin(), pvalues.end())));
            //     cout << "Excluding file: " << inFileListAux[std::distance(pvalues.begin(),std::min_element(pvalues.begin(), pvalues.end()))] << endl;
            //     inFileListAux.erase(inFileListAux.begin()+std::distance(pvalues.begin(),std::min_element(pvalues.begin(), pvalues.end())));
            //     if (inFileListAux.size()==1)
            //     {
            //         break;
            //     }
            // }
            // else
            // {
            //     go=false;
            //     h1->Draw();
            // }
            
            

            


            int position = std::distance(pvalues.begin(),std::min_element(pvalues.begin(), pvalues.end()));
            // one-tail exclusion
            cout << "------------" << endl;
            cout << pvalues[position] << endl;
            cout << "------------" << endl;
            // if (pvalues[std::distance(pvalues.begin(),std::min_element(pvalues.begin(), pvalues.end()))]>0.0001)
            if (pvalues[position]>0.3/pvalues.size())
            {
                go=false;
                h1->Draw();
            }
            else
            {
                myfile << inFileListAux[position]+","+std::to_string(rates[hdu-1][position])+","+std::to_string(ratesErr[hdu-1][position])+","+std::to_string(pvalues[position])+","+std::to_string(two_e_1pix[hdu-1][position])+","+std::to_string(two_e_h[hdu-1][position])+","+std::to_string(two_e_v[hdu-1][position])+","+std::to_string(two_e_d[hdu-1][position])+","+std::to_string(two_e_v[hdu-1][position]+two_e_d[hdu-1][position])+","+std::to_string(three_e[hdu-1][position])+",\n";
                cout << pvalues[position] << endl;
                rates[hdu-1].erase(rates[hdu-1].begin()+position);
                ratesErr[hdu-1].erase(ratesErr[hdu-1].begin()+position);
                entries[hdu-1].erase(entries[hdu-1].begin()+position);
                cout << "Excluding file: " << inFileListAux[position] << endl;
                inFileListAux.erase(inFileListAux.begin()+position);
                h1->Draw();
                canvas2->cd();
                canvas2->Print((outfilename+"_test.pdf").c_str());
                canvas2->Clear();
            }  


        }

        for (size_t i = 0; i < pvalues.size(); i++)
        {
            myfile << inFileListAux[i]+","+std::to_string(rates[hdu-1][i])+","+std::to_string(ratesErr[hdu-1][i])+","+std::to_string(pvalues[i])+","+std::to_string(two_e_1pix[hdu-1][i])+","+std::to_string(two_e_h[hdu-1][i])+","+std::to_string(two_e_v[hdu-1][i])+","+std::to_string(two_e_d[hdu-1][i])+","+std::to_string(two_e_v[hdu-1][i]+two_e_d[hdu-1][i])+","+std::to_string(three_e[hdu-1][i])+",\n";        
        }
        myfile << "*,*,*,\n";

        //// create experimental ladder
        sort(pvalues.begin(), pvalues.end());
        double ladderexp[int(pvalues.size())][3];
        for (size_t j = 0; j < pvalues.size(); j++)
        {
            ladderexp[j][0]=pvalues[j];
            ladderexp[j][1]=j/float(pvalues.size());
            ladderexp[j][2]=0; // 0 if data 1 if theory
        }
        


        // compute largest diff for experimental data
        double escalera[nsim+pvalues.size()][3];
        concatenate_doubles(escalera, ladder,ladderexp,nsim,pvalues.size());
        qsort(escalera, ARRSIZE(escalera), sizeof(*escalera), compare);
        double largest_diff=0;
        compute_largest_diff(escalera, pvalues, largest_diff);
        double real_largest_diff=largest_diff;

        // compute distribution of largest difference to calculate KS p-value
        double counter=0;
        size_t npvalue=10000;
        auto rndm = TRandom3(0);
	    rndm.SetSeed(0);
        for (size_t i = 0; i < npvalue; i++) // our p-value will be exact as 1/10000
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

            if (largest_diff>real_largest_diff)
            {
                counter++;
            }
              
        }

        char tv_KS[40];
        sprintf(tv_KS,  "p-value = %3.2f",counter/npvalue);
        TPaveLabel *px = new TPaveLabel(0.62,0.15,0.83,0.24,tv_KS,"brNDC");
        px->SetBorderSize(0);px->SetTextFont(42); px->SetTextSize(0.5);px->SetTextAlign(22);px->SetTextColor(1); 
        px->Draw();
        cout << "p-value: " << counter/npvalue << endl;
        
        



    }
    
    canvas2->cd();
    canvas2->Print((outfilename+"_test.pdf").c_str());
    canvas2->Clear();
    canvas2->Print((outfilename+"_test.pdf]").c_str());

    TCanvas * canvas3   = new TCanvas("Plots","Plots",0,0,1200,1000); //coord top lef, size x, sizey
    canvas3->Print((outfilename+"_rates.pdf[").c_str());
    // canvas3->Divide(2,2);
    outfile->cd();
    // TF1 fitfunc2("f","[0]*(TMath::Gaus(x,[1],[2],1))",0,1000);
    // TF1 * fitfunc2 = new TF1("f","[0]*(TMath::Gaus(x,[1],[2],1))",0,1000);
    vector <TF1*> f0(4);
    for (size_t i = 0; i < 4; i++)
    {
        TF1 *a = new TF1(("f"+std::to_string(i)).c_str(),"[0]*(TMath::Gaus(x,[1],[2],1))",0,1);
        f0[i] = (a);
    }

    //range
    float pixmax = 5e-3;
    float pixmin = 0;
    int pixnbin = 100;
    // double pixbinning = "({0},{1},{2})".format(pixnbin,pixmin,pixmax); std::to_string(pixnbin)+","+
    float binwidth = (pixmax-pixmin)/pixnbin;
    gStyle->SetOptStat(1);
    for (size_t hdu = 1; hdu < 5; hdu++)
    {

        TH1F * h3 = new TH1F("h3","h3",pixnbin,pixmin,pixmax);
        // canvas3->cd(hdu);
        for (size_t i = 0; i < rates[hdu-1].size(); i++)
        {
            h3->Fill(rates[hdu-1][i]);
        }
        h3->Draw();

        f0[hdu-1]->SetParameters(h3->GetEntries()*binwidth,h3->GetMean(),h3->GetStdDev());
        h3->Fit(f0[hdu-1],"QSL","",pixmin,pixmax);
        // TF1 *fit = h3->GetFunction(f0[hdu-1]->GetName());
        TF1 *fit = h3->GetFunction(("f"+std::to_string(hdu-1)).c_str());
        // TF1 *fit = h0->GetFunction("fitfunc");
        f0[hdu-1]->SetParameters(fit->GetParameter(0),fit->GetParameter(1),fit->GetParameter(2));
        cout << fit->GetParameter(0) << endl;
        cout << f0[hdu-1]->GetParameter(0) << endl;
        cout << fit->GetParameter(1) << endl;
        cout << f0[hdu-1]->GetParameter(1) << endl;
        cout << fit->GetParameter(2) << endl;
        cout << f0[hdu-1]->GetParameter(2) << endl;
        f0[hdu-1]->SetNpx(100000);
        canvas3->cd();
        canvas3->Clear();

    
    }

    auto legend = new TLegend(0.60,0.25,0.85,0.40); 
    canvas3->cd();
    for (size_t hdu = 1; hdu < 5; hdu++)
    {
        f0[hdu-1]->SetLineColor(hdu);
        legend->AddEntry(f0[hdu-1],("Hdu"+std::to_string(hdu)).c_str(),"f");
        if (hdu==1)
        {
            f0[hdu-1]->Draw("");
            f0[hdu-1]->GetYaxis()->SetRangeUser(0,15);
            f0[hdu-1]->GetXaxis()->SetRangeUser(pixmin,pixmax);
            f0[hdu-1]->GetXaxis()->SetTitle("mu");
            f0[hdu-1]->SetTitle("1e- events");
            f0[hdu-1]->Draw("");
        }
        else
        {
            f0[hdu-1]->Draw("SAME");
        }
        legend->Draw();
    }

    canvas3->cd();
    canvas3->Print((outfilename+"_rates.pdf").c_str());
    canvas3->Clear();
    canvas3->Print((outfilename+"_rates.pdf]").c_str());
    
    outfile->Write();
    outfile->Close();
    myfile.close();

    
    return 0;
}