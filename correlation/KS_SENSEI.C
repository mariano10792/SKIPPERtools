/*
Sanity Checks for SENSEI Data processed by Mariano Cababie
Date: October 1, 2020
Author: Dario Rodrigues
email: dariorodriguesfm@gmail.com

*/

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include "TObject.h"
#include "TKey.h"
#include "TString.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include <numeric>
#include "fitsio.h"
#include <ctime>
#include <cstdlib>
#include <chrono>
#include <sstream>
#include "TH1D.h"
#include "TH2D.h"
#include "TStopwatch.h"
#include <TROOT.h>
#include <TStyle.h>
#include "TMath.h"
#include "TGraphErrors.h"
#include "TLeafC.h"
#include "TStyle.h"
#include "TChain.h"
#include "TPaveLabel.h"
#include "TLegend.h"
#include <iterator>
#include <algorithm>
// #include "Minuit2/Minuit2Minimizer.h"
#include "Math/Functor.h"
#include "TCanvas.h"
#include "TPad.h"
#include <stdio.h>
#include <stdlib.h>

#define ARRSIZE(arr) (sizeof(arr)/sizeof(*(arr)))

using namespace std;

template <class Container>

void split(const std::string& str, Container& cont, char delim = ' ')
{
    std::stringstream ss(str);
    std::string token;
    while (std::getline(ss, token, delim)) {
        cont.push_back(token);
    }
}



//// Selection options ////

////////////////////////////////////////////////////////////////////////////////////

int runID1e_prev;
int N;
int x; 
int y; 
int ohdu; 
int runID;
int nsim;
int factor;

vector<vector<double>> pvalores(4);


//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////


void Enable_and_Set_Branches_Experimental(TTree* & tree){

	tree->SetBranchStatus("*",1); //enable all branches

	
	tree->SetBranchAddress("RUNID",&runID);
	tree->SetBranchAddress("ohdu",&ohdu);
	tree->SetBranchAddress("x",&x);
	tree->SetBranchAddress("y",&y);
	
}


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


//This function sorts values and creates x axis for each value
void sorting(vector<vector<double>> &xaxis, vector<vector<double>> &DATA)
{
	int cExp = DATA[3].size();
	int nExp = DATA[0].size();

	sort(DATA[0].begin(), DATA[0].end());
	sort(DATA[1].begin(), DATA[1].end());
	sort(DATA[2].begin(), DATA[2].end());
	sort(DATA[3].begin(), DATA[3].end());

	for (Int_t i=0;i<nExp;i++) {
		xaxis[0].push_back((i+1)*1./nExp); 	
	}

	for (Int_t i=0;i<cExp;i++) {
		xaxis[1].push_back((i+1)*1./cExp);
	}

}


void analyse(vector<vector<double>> &events, int Entries, TTree * tree, int hdu)
{
	Enable_and_Set_Branches_Experimental(tree);
	vector<int> runID1e;  
	vector<int> ohdu1e;

	for(int i_event=0;i_event<Entries; i_event++){

		tree->GetEntry(i_event);
		if(ohdu!=hdu){continue;}
				
		events[0].push_back(x);
		events[1].push_back(y);
		events[2].push_back(sqrt(pow(x,2)+pow(y,2)));

		runID1e.push_back(runID);
		ohdu1e.push_back(ohdu);
	}

	const int n = events[0].size();

	double dx, dy, distance;
	int neighbours=3;

	vector <float_t> distancetemp(neighbours,3200);
	for (int i = 0; i < n; ++i){
		distancetemp.assign(distancetemp.size(),3200);
		if (i>=n-neighbours){distancetemp.pop_back();}
		for (int j = i+1; j < n; ++j){
			if (runID1e[i]==runID1e[j]){
				if (ohdu1e[i]==ohdu1e[j]){
					dx=events[0][i]-events[0][j];
					dy=events[1][i]-events[1][j];
					distance=sqrt(pow(dx,2)+pow(dy,2));
					if (distance<distancetemp[distancetemp.size()-1])
					{
						distancetemp[distancetemp.size()-1]=distance;
						sort(distancetemp.begin(), distancetemp.end());
					}
				}
			}
		}
		for (size_t temp = 0; temp < distancetemp.size(); temp++)
		{
			events[3].push_back(distancetemp[temp]);
		}
		events[3].pop_back();
	}
	const int c = events[3].size();
}


void kolmogorovTest(vector<vector<double>> ON, vector<vector<double>> OFF, vector<vector<double>> x_ON, vector<vector<double>> x_OFF, int nn, int mm, double &pvalue, int index1, int index2)
{
	int ntotal=ON[0].size()+OFF[0].size();
    int ctotal=ON[3].size()+OFF[3].size(); //4th test

    double x_ONOFF[ntotal][3];
    double y_ONOFF[ntotal][3];
    double d_ONOFF[ntotal][3];
    double c_ONOFF[ctotal][3]; //4th test
	double diff,last_ON,last_OFF,largest_diff=0;
	for (Int_t i=0;i<nn;i++) {
	x_ONOFF[i][0]=ON[index1][i];
	x_ONOFF[i][1]=x_ON[index2][i];
	x_ONOFF[i][2]=1;
	}
	for (Int_t i=0;i<mm;i++) {
		x_ONOFF[i+nn][0]=OFF[index1][i];
		x_ONOFF[i+nn][1]=x_OFF[index2][i];
		x_ONOFF[i+nn][2]=0;
	}
	qsort(x_ONOFF, ARRSIZE(x_ONOFF), sizeof(*x_ONOFF), compare);

	for (Int_t i=0;i<ntotal;i++) {
		if (x_ONOFF[i][2]==0) last_OFF = x_ONOFF[i][1];
		if (x_ONOFF[i][2]==1) last_ON = x_ONOFF[i][1];
		diff=abs(last_ON-last_OFF);
		if (largest_diff<diff) largest_diff=diff;
	}

	// cout<<endl;
	// cout<<"Kolmogorov–Smirnov statistic for X-coordinate is: "<<largest_diff<<endl;

	// From https://en.wikipedia.org/wiki/Kolmogorov%E2%80%93Smirnov_test
	// double pvalue_x_old=2*exp(-2*largest_diff*largest_diff*(nn*mm)/(nn+mm));
	// cout<<"p-value_old X: "<<pvalue_x_old<<endl;

	pvalue=0;
	double sum_x=0;
	double z_x=largest_diff*pow((nn*mm)*1./(nn+mm),0.5);
	for (int r = 1; r < 100; ++r) pvalue+=2*pow(-1,r-1)*exp(-2*pow(r,2)*pow(z_x,2));

	pvalores[index1].push_back(pvalue);
}

int kolmogorov(vector<vector<double>> OFF, vector<vector<double>> ON, int hdu, std::string outfilename){
	

//////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// Actually Kolmogorov starts here ////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
	
	vector<vector<double>> x_ON(2);
	sorting(x_ON,ON);
	vector<vector<double>> x_OFF(2);
	sorting(x_OFF,OFF);

//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////// KOLMOGOROV //////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////

	int nn=ON[0].size();
	int mm=OFF[0].size();   

    // double diff,last_ON,last_OFF,largest_diff;

	double pvalue_x=0; double pvalue_y=0; double pvalue_c=0; double pvalue_d=0;

	kolmogorovTest(ON, OFF, x_ON, x_OFF, nn, mm, pvalue_x,0,0);
	kolmogorovTest(ON, OFF, x_ON, x_OFF, nn, mm, pvalue_y,1,0);
	kolmogorovTest(ON, OFF, x_ON, x_OFF, nn, mm, pvalue_d,2,0);
	kolmogorovTest(ON, OFF, x_ON, x_OFF, nn, mm, pvalue_c,3,1);


//////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////// PLOTS ////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
	bool PLOTS=true;

	if (PLOTS){

	TCanvas * canvas   = new TCanvas("Plots","Plots",0,0,2000,1000); //coord top lef, size x, sizey
	canvas->Divide(2,2);

	canvas->cd(1);
	canvas->cd(1)->SetGridx();
	canvas->cd(1)->SetGridy();

	// TGraph * grxxON = new TGraph(nExp_ON,xxx_ON,x_ON);
	TGraph * grxxON = new TGraph(ON[0].size(), &ON[0][0], &x_ON[0][0]);
    grxxON->Draw("");
	grxxON->GetHistogram()->GetXaxis()->SetLimits(0,3200);
	grxxON->Draw("");
    grxxON->GetXaxis()->SetTitle("X-coordinate");
    grxxON->GetYaxis()->SetTitle("Cumulative probability");
	
	// grxxON->SetMarkerStyle(21);
	grxxON->SetMarkerColor(kBlue);
	grxxON->SetLineColor(kBlue);
	grxxON->SetLineWidth(2);
	grxxON->SetTitle("Empirical distribution functions");	

////////////////////////////////////////////////////////////////

	// TGraph * grxxOFF = new TGraph(nExp_OFF,xxx_OFF,x_OFF);
	TGraph * grxxOFF = new TGraph(OFF[0].size(), &OFF[0][0], &x_OFF[0][0]);
    grxxOFF->SetLineColor(kRed);
    grxxOFF->SetLineWidth(2);
    grxxOFF->Draw("same");
    
	// Set limit y axis
	// gryy->GetHistogram()->SetMinimum(-0.1);
	// gryy->GetHistogram()->SetMaximum(1.1); 

	// Set limit x axis

	auto legend = new TLegend(0.60,0.25,0.85,0.40); 
	legend->AddEntry(grxxON,"Real Data","f");
	legend->AddEntry(grxxOFF,"Synthetic Data","f");
	legend->Draw();

    char tv_KS[40];
    sprintf(tv_KS,  "p-value = %3.2f",pvalue_x);
    TPaveLabel *px = new TPaveLabel(0.62,0.15,0.83,0.24,tv_KS,"brNDC");
    px->SetBorderSize(0);px->SetTextFont(42); px->SetTextSize(0.5);px->SetTextAlign(22);px->SetTextColor(1); 
    px->Draw();


//////////////////////////////////////////////////////////////////////////////////////////////////////

	canvas->cd(2);
	canvas->cd(2)->SetGridx();
	canvas->cd(2)->SetGridy();

	// TGraph * gryyON = new TGraph(nExp_ON,yyy_ON,y_ON);
	TGraph * gryyON = new TGraph(ON[1].size(), &ON[1][0], &x_ON[0][0]);
	
    gryyON->Draw("");
	gryyON->GetHistogram()->GetXaxis()->SetLimits(0,500);
	gryyON->Draw("");
    gryyON->GetXaxis()->SetTitle("Y-coordinate");
    gryyON->GetYaxis()->SetTitle("Cumulative probability");
	
	// gryyON->SetMarkerStyle(21);
	gryyON->SetMarkerColor(kBlue);
	gryyON->SetLineColor(kBlue);
	gryyON->SetLineWidth(2);
	gryyON->SetTitle("Empirical distribution functions");	

////////////////////////////////////////////////////////////////

	// TGraph * gryyOFF = new TGraph(nExp_OFF,yyy_OFF,y_OFF);
	TGraph * gryyOFF = new TGraph(OFF[1].size(), &OFF[1][0], &x_OFF[0][0]);
    gryyOFF->SetLineColor(kRed);
    gryyOFF->SetLineWidth(2);
	// gryyOFF->SetAxisRange(0., 500.,"x");
    gryyOFF->Draw("same");
    
	// Set limit y axis
	// gryy->GetHistogram()->SetMinimum(-0.1);
	// gryy->GetHistogram()->SetMaximum(1.1); 

	// Set limit x axis

	legend->Draw();

	// char tv_KS[40];
    sprintf(tv_KS,  "p-value = %3.2f",pvalue_y);
    TPaveLabel *py = new TPaveLabel(0.62,0.15,0.83,0.24,tv_KS,"brNDC");
    py->SetBorderSize(0);py->SetTextFont(42); py->SetTextSize(0.5);py->SetTextAlign(22);py->SetTextColor(1); 
    py->Draw();

//////////////////////////////////////////////////////////////////////////////////////////////////////

	canvas->cd(3);
	canvas->cd(3)->SetGridx();
	canvas->cd(3)->SetGridy();

	// TGraph * grddON = new TGraph(nExp_ON,ddd_ON,d_ON);
	TGraph * grddON = new TGraph(ON[2].size(), &ON[2][0], &x_ON[0][0]);
	
	
    grddON->Draw("");
	grddON->GetHistogram()->GetXaxis()->SetLimits(0,3200);
	grddON->Draw("");
    grddON->GetXaxis()->SetTitle("Distance to origin");
    grddON->GetYaxis()->SetTitle("Cumulative probability");
	
	// grddON->SetMarkerStyle(21);
	grddON->SetMarkerColor(kBlue);
	grddON->SetLineColor(kBlue);
	grddON->SetLineWidth(2);
	grddON->SetTitle("Empirical distribution functions");

////////////////////////////////////////////////////////////////

	// TGraph * grddOFF = new TGraph(nExp_OFF,ddd_OFF,d_OFF);
	TGraph * grddOFF = new TGraph(OFF[2].size(), &OFF[2][0], &x_OFF[0][0]);
    grddOFF->SetLineColor(kRed);
    grddOFF->SetLineWidth(2);
	grddOFF->GetXaxis()->SetRange(0,3200);
    grddOFF->Draw("same");
    
	// Set limit y axis
	// gryy->GetHistogram()->SetMinimum(-0.1);
	// gryy->GetHistogram()->SetMaximum(1.1); 

	// Set limit x axis

	legend->Draw();

	// char tv_KS[40];
    sprintf(tv_KS,  "p-value = %3.2f",pvalue_d);
    TPaveLabel *pd = new TPaveLabel(0.62,0.15,0.83,0.24,tv_KS,"brNDC");
    pd->SetBorderSize(0);pd->SetTextFont(42); pd->SetTextSize(0.5);pd->SetTextAlign(22);pd->SetTextColor(1); 
    pd->Draw();

//////////////////////////////////////////////////////////////////////////////////////////////////////

	canvas->cd(4);
	canvas->cd(4)->SetGridx();
	canvas->cd(4)->SetGridy();

	// TGraph * grccON = new TGraph(cExp_ON,ccc_ON,c_ON);
	TGraph * grccON = new TGraph(ON[3].size(), &ON[3][0], &x_ON[1][0]);
	
    grccON->Draw("");
	grccON->GetHistogram()->GetXaxis()->SetLimits(0,3200);
	grccON->Draw("");
    grccON->GetXaxis()->SetTitle("Distance between events");
    grccON->GetYaxis()->SetTitle("Cumulative probability");
	
	// grccON->SetMarkerStyle(21);
	grccON->SetMarkerColor(kBlue);
	grccON->SetLineColor(kBlue);
	grccON->SetLineWidth(2);
	grccON->SetTitle("Empirical distribution functions");	

////////////////////////////////////////////////////////////////

	// TGraph * grccOFF = new TGraph(cExp_OFF,ccc_OFF,c_OFF);
	TGraph * grccOFF = new TGraph(OFF[3].size(), &OFF[3][0], &x_OFF[1][0]);
    grccOFF->SetLineColor(kRed);
    grccOFF->SetLineWidth(2);
    grccOFF->Draw("same");
	gPad->SetLogx(1);
    
	// Set limit y axis
	// gryy->GetHistogram()->SetMinimum(-0.1);
	// gryy->GetHistogram()->SetMaximum(1.1); 

	// Set limit x axis

	legend->Draw();

	// char tv_KS[40];
    sprintf(tv_KS,  "p-value = %3.2f",pvalue_c);
    TPaveLabel *pc = new TPaveLabel(0.62,0.15,0.83,0.24,tv_KS,"brNDC");
    pc->SetBorderSize(0);pc->SetTextFont(42); pc->SetTextSize(0.5);pc->SetTextAlign(22);pc->SetTextColor(1); 
    pc->Draw();
	canvas->SaveAs(("files/KS_"+outfilename+"_hdu"+std::to_string(hdu)+".pdf").c_str());
	}
	return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[]){

	gROOT->SetStyle("Plain");
	cout.precision(4);  // Print using four decimals precision
	cout<<fixed;        // ...exactly four decimals
	
	//set run to analyse (argument)
	size_t initrun = atoi(argv[argc-1]);

	//get inFileList
	std::vector<string> inFileList;
	for(int i=optind; i<argc-1; ++i){ // -1 to delete option -h or -p. Not final solution.
		inFileList.push_back(argv[i]);
	}

	//outfilename
	std::vector<std::string> words;
	split(inFileList[0],words,'/');
	std::string outfilename;
	outfilename=words.back();
	outfilename="pvalues_run"+std::to_string(initrun)+"_"+outfilename;

	//initialize output tree
	TFile *outputFile = new TFile(("files/"+outfilename).c_str(),"RECREATE"); //output file
	cout << outputFile->GetName() << endl;
	Float_t pvalues1Val; Float_t pvalues2Val; Float_t pvalues3Val; Float_t pvalues4Val; Int_t ohduVal; Int_t runIDVal; Int_t LTANAMEVal;
	TTree * pvaluesTree = new TTree("pvalues","pvalues");
	pvaluesTree->Branch("pvalues1",&pvalues1Val,"pvalues1/F");
	pvaluesTree->Branch("pvalues2",&pvalues2Val,"pvalues2/F");
	pvaluesTree->Branch("pvalues3",&pvalues3Val,"pvalues3/F");
	pvaluesTree->Branch("pvalues4",&pvalues4Val,"pvalues4/F");
	pvaluesTree->Branch("runID",&runIDVal,"runID/I");
	pvaluesTree->Branch("ohdu",&ohduVal,"ohdu/I");
	pvaluesTree->Branch("LTANAME",&LTANAMEVal,"LTANAME/I");


	for (size_t iFile = 0; iFile < inFileList.size(); iFile++)
	{
		
		// Open Experimental file and tree
		TFile * f_ON = TFile::Open(inFileList[iFile].c_str());
		cout << "Opening " << inFileList[iFile] << endl;
		if (!f_ON->IsOpen()) {std::cerr << "ERROR: cannot open the root file with experimental data" << std::endl;}
		TTree * t_ON=nullptr;
		t_ON = (TTree*) f_ON->Get("calPixTree1e");
		int Entries_ON = t_ON -> GetEntries();
		cout<<"Entries in data file from Experimental images: "<<Entries_ON<<endl;

		// Open Simulated file
		TFile * f_OFF = TFile::Open(inFileList[iFile].c_str());
		if (!f_OFF->IsOpen()) {std::cerr << "ERROR: cannot open the root file with experimental data" << std::endl;}

		// Open Simulated header Tree
		TTree * simHeader = (TTree*) f_OFF->Get("simHeader");
		int nsim; int factor; int LTANAME;
		simHeader->SetBranchStatus("*",1);
		simHeader->SetBranchAddress("factor",&factor);
		simHeader->SetBranchAddress("nsim",&nsim);
		simHeader->SetBranchAddress("LTANAME",&LTANAME);
		simHeader->GetEntry(0);
		cout << "Number of simulations: " << nsim << endl;
		cout << "Factor of events simulated: " << factor<< endl;

		




		//analyse and compare data
		for (int hdu = 1; hdu < 5; hdu++){


			//Get number of entries for each run hdu

			vector <int> entries;
			int Entries_ON = t_ON -> GetEntries(("ohdu=="+std::to_string(hdu)).c_str());
			entries.push_back(Entries_ON);
			for (size_t run = initrun; run < initrun+1; run++)
			{
				TTree * t_OFF=nullptr;
				t_OFF = (TTree*) f_OFF->Get(("simPixTree"+std::to_string(run)).c_str());
				int Entries_OFF = t_OFF -> GetEntries(("ohdu=="+std::to_string(hdu)).c_str());
				entries.push_back(Entries_OFF);
			}
			
			int minEntries = *std::min_element(entries.begin(), entries.end());
			cout << minEntries << endl;




			// analyse Experimental data 
			vector<vector<double>> ON(4);
			Entries_ON = t_ON->GetEntries();
			analyse(ON, Entries_ON, t_ON, hdu);

			// analyse Simulated data and perform Kolmogorov test
			for (size_t run = initrun; run < initrun+1; run++)
			{
				vector<vector<double>> OFF(4);
				TTree * t_OFF=nullptr;
				t_OFF = (TTree*) f_OFF->Get(("simPixTree"+std::to_string(run)).c_str());
				int Entries_OFF = t_OFF -> GetEntries();
				cout<<"Entries in data file from Simulated images: " << Entries_OFF <<endl;
				analyse(OFF, Entries_OFF, t_OFF, hdu);
				kolmogorov(OFF, ON, hdu, outfilename);
			}



			
			//write in output tree
			ohduVal=hdu;
			cout << ohduVal << endl;
			runIDVal=runID;
			cout << runIDVal << endl;
			LTANAMEVal=LTANAME;
			cout << LTANAMEVal << endl;
			pvalues1Val= std::accumulate( pvalores[0].begin(), pvalores[0].end(), 0.0) / pvalores[0].size(); 
			pvalues2Val= std::accumulate( pvalores[1].begin(), pvalores[1].end(), 0.0) / pvalores[1].size();
			pvalues3Val= std::accumulate( pvalores[2].begin(), pvalores[2].end(), 0.0) / pvalores[2].size();
			pvalues4Val= std::accumulate( pvalores[3].begin(), pvalores[3].end(), 0.0) / pvalores[3].size();
			cout << pvalues1Val << endl;
			cout << pvalues2Val << endl;
			cout << pvalues3Val << endl;
			cout << pvalues4Val << endl;
			outputFile->cd();
			pvaluesTree->Fill();
			  
			pvalores[0].clear(); pvalores[1].clear(); pvalores[2].clear(); pvalores[3].clear();
			
		}

	}
	pvaluesTree->Write();
	outputFile->Write();
	outputFile->Close();
	return 0;
}


