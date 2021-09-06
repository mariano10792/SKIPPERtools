/*
Sanity Checks for SENSEI Data processed by Mariano Cababie
Date: October 1, 2020
Author: Dario Rodrigues
email: dariorodriguesfm@gmail.com

*/

#include <fstream>
#include <iostream>
#include <random>
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


// std::default_random_engine generator;
// std::uniform_real_distribution<double> distribution(0.0,116000.0);

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
vector<vector<double>> pvaloresbis(4);
vector<vector<double>> largest_diff_sim(4);


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


// This function counts pvalues
void count(Float_t &largestdiffVal, vector<vector<double>> largest_diff_sim, int index1)
{	
	largestdiffVal=0;
	cout << largest_diff_sim[index1].size() << endl;
	cout << "[";
	cout << largest_diff_sim[index1][0] << ",";
	for (size_t i = 1; i < largest_diff_sim[index1].size(); i++)
	{
		if (largest_diff_sim[index1][0]>largest_diff_sim[index1][i])
		{
			largestdiffVal++;
		}
		if (largest_diff_sim[index1][0]==largest_diff_sim[index1][i])
		{
			largestdiffVal=largestdiffVal+0.5;
		}
		
		cout << largest_diff_sim[index1][i] <<",";
	}
	cout << "]" << endl;
	largestdiffVal=largestdiffVal/largest_diff_sim[index1].size();
	cout << largestdiffVal << endl;
}


//This function sorts values and creates x axis for each value
void sorting(vector<vector<double>> &xaxis, vector<vector<double>> &DATA)
{
	int nExp = DATA[0].size(); //distnace from x origin
	int cExp = DATA[3].size(); //distance entre pares

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


void analyse(vector<vector<double>> &events, int Entries, TTree * tree, int hdu, int neighbours)
{
	Enable_and_Set_Branches_Experimental(tree);
	vector<int> runID1e;  
	vector<int> ohdu1e;
	// auto test = TRandom3(0);
	// test.SetSeed(0);
	for(int i_event=0;i_event<Entries; i_event++){

		tree->GetEntry(i_event);
		if(ohdu!=hdu){continue;}
		
		events[0].push_back(x);
		events[1].push_back(y);
		events[2].push_back(sqrt(pow(x,2)+pow(y,2)));

		runID1e.push_back(runID);
		ohdu1e.push_back(ohdu);
	}

	// for(int i_event=0;i_event<Entries; i_event++){

	// 	tree->GetEntry(i_event);
	// 	if(ohdu!=hdu){continue;}
	// 	double  xx= test.Uniform(0,1)*1160.0;
	// 	events[0].push_back(xx);
	// 	runID1e.push_back(runID);
	// 	ohdu1e.push_back(ohdu);
	// }
	// test.SetSeed(0);
	// for(int i_event=0;i_event<Entries; i_event++){
	// 	if(ohdu!=hdu){continue;}
	// 	double yy= test.Uniform(0,1)*1160.0;
	// 	events[1].push_back(yy);
	// }
	// for(int i_event=0;i_event<Entries; i_event++){
	// 	events[2].push_back(sqrt(pow(events[0][i_event],2)+pow(events[1][i_event],2)));
	// }

	const int n = events[0].size();

	// double dx, dy, distance;

	// vector <float_t> distancetemp(neighbours,3200);
	// for (int i = 0; i < n; ++i){
	// 	distancetemp.assign(distancetemp.size(),3200);
	// 	// for (size_t j = 0; j < distancetemp.size(); j++)
	// 	// {
	// 	// 	cout << distancetemp[j] << endl;
	// 	// }
	// 	// cout << ".................." << endl;

	// 	if (i>=n-neighbours){distancetemp.pop_back();}
	// 	for (int j = i+1; j < n; ++j){
	// 		if (runID1e[i]==runID1e[j]){
	// 			if (ohdu1e[i]==ohdu1e[j]){
	// 				dx=events[0][i]-events[0][j];
	// 				dy=events[1][i]-events[1][j];
	// 				distance=sqrt(pow(dx,2)+pow(dy,2));
	// 				if (distance<distancetemp[distancetemp.size()-1])
	// 				{
	// 					distancetemp[distancetemp.size()-1]=distance;
	// 					// cout << distance << endl;
	// 					sort(distancetemp.begin(), distancetemp.end());
	// 				}
	// 			}
	// 		}
	// 	}

	// 	// for (size_t j = 0; j < distancetemp.size(); j++)
	// 	// {
	// 	// 	cout << distancetemp[j] << endl;
	// 	// }
	// 	// cout << "------------------" << endl;
	// 	for (size_t temp = 0; temp < distancetemp.size(); temp++)
	// 	{
	// 		events[3].push_back(distancetemp[temp]);
	// 	}
	// 	events[3].pop_back();
	// }

	// cout << "*****************" << endl;
	// cout << "*****************" << endl;
	// cout << "*****************" << endl;
	// cout << "*****************" << endl;
	
	const int c = events[3].size();
	cout << c << " pair of events counted." << endl;
}


void kolmogorovTest2(vector<vector<double>> ON, vector<vector<double>> OFF, vector<vector<double>> x_ON, vector<vector<double>> x_OFF, double &pvalue, int index1, int index2)
{

	// int nn=10000;
	// int mm=OFF[0].size();
	// if (index2==1){
	// 	nn=100000; mm=OFF[3].size();
	// }
		
	// int ntotal=nn+mm;
	
    // double x_ONOFF[ntotal][3];

	// double diff,last_ON,last_OFF,largest_diff=0;
	// for (Int_t i=0;i<nn;i++) {
	// 	if (index1==0) x_ONOFF[i][0]=(double(i)/nn)*1160;
	// 	if (index1==1) x_ONOFF[i][0]=(double(i)/nn)*1160;
	// 	x_ONOFF[i][1]=double(i)/nn;
	// 	x_ONOFF[i][2]=1;
	// }


	int nn=ON[0].size();
	int mm=OFF[0].size();
	if (index2==1){
		nn=ON[3].size(); mm=OFF[3].size();
	}
	int ntotal=nn+mm;

    double x_ONOFF[ntotal][3];


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
	largest_diff_sim[index1].push_back(largest_diff);
	

}

void kolmogorovTest(vector<vector<double>> ON, vector<vector<double>> OFF, vector<vector<double>> x_ON, vector<vector<double>> x_OFF, double &pvalue, int index1, int index2)
{
	int nn=ON[0].size();
	int mm=OFF[0].size();
	if (index2==1){
		nn=ON[3].size(); mm=OFF[3].size();
	}
	int ntotal=nn+mm;
	
    // int ctotal=ON[3].size()+OFF[3].size(); //4th test

    double x_ONOFF[ntotal][3];
    // double y_ONOFF[ntotal][3];
    // double d_ONOFF[ntotal][3];
    // double c_ONOFF[ctotal][3]; //4th test

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
		// if (index2==0) cout << *x_ONOFF[i] <<endl;
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

	pvalue=0; double pvalue_total=0;
	double sum_x=0;
	// double z_x=largest_diff*pow((nn*mm)*1./(nn+mm),0.5);
	double z_x=exp(log(largest_diff)+0.5*(log(nn)+log(mm)-log(nn+mm)));
	// z_x = z_x+1/(6*pow(nn,0.5))+(z_x-1)/(4*(nn));
	// int r=1;
	// while (pvalue==0 || abs(pvalue)>0.0001)
	// {
	// 	pvalue=2*pow(-1,r-1)*exp(-2*pow(r,2)*pow(z_x,2));
	// 	cout << pvalue << endl;
	// 	pvalue_total+=pvalue;
	// 	r++;
	// 	if (r>10000)
	// 	{
	// 		cout << "FATAL ERROR: P-VALUE DOES NOT CONVERGE." << endl;
	// 		break;
	// 	}
		
	// }

	// cout << "It took me " << r << " terms to converge to a p-value." << endl; 
	// cout << "---------------------" << endl;
	

	for (int r = 1; r < 100; ++r) pvalue_total+=2*pow(-1,r-1)*exp(-2*pow(r,2)*pow(z_x,2));

	if (pvalue_total<0.5)
	{
		pvalores[index1].push_back(1-pvalue_total);
	}
	else
	{
		pvalores[index1].push_back(pvalue_total);
	}
	
	
	// pvalores[index1].push_back(pvalue_total);
	pvaloresbis[index1].push_back(pvalue_total);
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

	// kolmogorovTest(ON, OFF, x_ON, x_OFF, pvalue_x,0,0);
	kolmogorovTest2(ON, OFF, x_ON, x_OFF, pvalue_x,0,0);
	// kolmogorovTest(ON, OFF, x_ON, x_OFF, pvalue_y,1,0);
	kolmogorovTest2(ON, OFF, x_ON, x_OFF, pvalue_y,1,0);
	// kolmogorovTest(ON, OFF, x_ON, x_OFF, pvalue_d,2,0);
	kolmogorovTest2(ON, OFF, x_ON, x_OFF, pvalue_d,2,0);
	// kolmogorovTest(ON, OFF, x_ON, x_OFF, pvalue_c,3,1);
	kolmogorovTest2(ON, OFF, x_ON, x_OFF, pvalue_c,3,1);


//////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////// PLOTS ////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
	bool PLOTS=false;

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
	// grxxON->GetHistogram()->GetXaxis()->SetLimits(0,1160);
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
	// gryyON->GetHistogram()->GetXaxis()->SetLimits(0,1160);
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
	// grddON->GetHistogram()->GetXaxis()->SetLimits(0,1160);
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
	// grddOFF->GetXaxis()->SetRange(0,1160);
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
	// grccON->GetHistogram()->GetXaxis()->SetLimits(0,1160);
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

	TStopwatch t; // time counter
	t.Start();

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
	Int_t ohduVal; Int_t runIDVal; Int_t LTANAMEVal;
	// Float_t pvalues1Val; Float_t pvalues2Val; Float_t pvalues3Val; Float_t pvalues4Val; 
	// Float_t pvaluesbis1Val; Float_t pvaluesbis2Val; Float_t pvaluesbis3Val; Float_t pvaluesbis4Val;
	// TTree * pvaluesTree = new TTree("pvalues","pvalues");
	// pvaluesTree->Branch("pvalues1",&pvalues1Val,"pvalues1/F");
	// pvaluesTree->Branch("pvalues2",&pvalues2Val,"pvalues2/F");
	// pvaluesTree->Branch("pvalues3",&pvalues3Val,"pvalues3/F");
	// pvaluesTree->Branch("pvalues4",&pvalues4Val,"pvalues4/F");
	// pvaluesTree->Branch("pvaluesbis1",&pvaluesbis1Val,"pvaluesbis1/F");
	// pvaluesTree->Branch("pvaluesbis2",&pvaluesbis2Val,"pvaluesbis2/F");
	// pvaluesTree->Branch("pvaluesbis3",&pvaluesbis3Val,"pvaluesbis3/F");
	// pvaluesTree->Branch("pvaluesbis4",&pvaluesbis4Val,"pvaluesbis4/F");
	// pvaluesTree->Branch("runID",&runIDVal,"runID/I");
	// pvaluesTree->Branch("ohdu",&ohduVal,"ohdu/I");
	// pvaluesTree->Branch("LTANAME",&LTANAMEVal,"LTANAME/I");
	TTree * largestdiffsTree = new TTree("largestdiffs","largestdiffs");
	Float_t largestdiff1Val; Float_t largestdiff2Val; Float_t largestdiff3Val; Float_t largestdiff4Val;
	largestdiffsTree->Branch("largestdiff1",&largestdiff1Val,"largestdiff1/F");
	largestdiffsTree->Branch("largestdiff2",&largestdiff2Val,"largestdiff2/F");
	largestdiffsTree->Branch("largestdiff3",&largestdiff3Val,"largestdiff3/F");
	largestdiffsTree->Branch("largestdiff4",&largestdiff4Val,"largestdiff4/F");
	largestdiffsTree->Branch("runID",&runIDVal,"runID/I");
	largestdiffsTree->Branch("ohdu",&ohduVal,"ohdu/I");
	largestdiffsTree->Branch("LTANAME",&LTANAMEVal,"LTANAME/I");


	for (size_t iFile = 0; iFile < inFileList.size(); iFile++)
	{
		
		// Open file
		TFile * file = TFile::Open(inFileList[iFile].c_str());
		cout << "Opening " << inFileList[iFile] << endl;
		if (!file->IsOpen()) {std::cerr << "ERROR: cannot open the root file with experimental data" << std::endl;}

		// Open Simulated header Tree
		TTree * simHeader = (TTree*) file->Get("simHeader");
		int nsim; int factor; int LTANAME;
		simHeader->SetBranchStatus("*",1);
		simHeader->SetBranchAddress("factor",&factor);
		simHeader->SetBranchAddress("nsim",&nsim);
		simHeader->SetBranchAddress("LTANAME",&LTANAME);
		simHeader->GetEntry(0);
		cout << "Number of simulations: " << nsim << endl;
		cout << "Factor of events simulated: " << factor<< endl;

		//set number of neighbours 
		int neighbours=1;


		//analyse and compare data
		for (int hdu = 1; hdu < 2; hdu++){



			// analyse "Theoretical" data (factor 100)
			vector<vector<double>> THEO(4); //x, y , x**2+y**2, distance entre eventos
			TTree * t_THEO=nullptr;
			t_THEO = (TTree*) file->Get("simPixTree00");
			// t_THEO = (TTree*) file->Get("simPixTree0");
			int Entries_THEO = t_THEO -> GetEntries();
			analyse(THEO, Entries_THEO, t_THEO, hdu, neighbours);


			// analyse Experimental data 
			vector<vector<double>> ON(4); //x, y , x**2+y**2, distance entre eventos
			TTree * t_ON=nullptr;
			t_ON = (TTree*) file->Get("calPixTree1e");
			int Entries_ON = t_ON -> GetEntries();
			analyse(ON, Entries_ON, t_ON, hdu , neighbours);

			// perform Kolmogorov test between Experimental and "Theoretical" data
			kolmogorov(THEO, ON, hdu, outfilename);

			// analyse Simulated data and perform Kolmogorov test
			// for (size_t run = 0; run < nsim-1; run++)
			// for (size_t run = 0; run < 100-1; run++)
			for (size_t run = 1; run < 100-1; run++)
			{
				vector<vector<double>> OFF(4);
				TTree * t_OFF=nullptr;
				t_OFF = (TTree*) file->Get(("simPixTree"+std::to_string(run)).c_str());
				int Entries_OFF = t_OFF -> GetEntries();
				cout<<"Entries in data file from Simulated images: " << Entries_OFF <<endl;
				analyse(OFF, Entries_OFF, t_OFF, hdu , neighbours);
				kolmogorov(THEO, OFF, hdu, outfilename);




			
				//write in output tree
				ohduVal=hdu;
				cout << ohduVal << endl;
				runIDVal=runID;
				cout << runIDVal << endl;
				LTANAMEVal=LTANAME;
				cout << LTANAMEVal << endl;
				// pvalues1Val= std::accumulate( pvalores[0].begin(), pvalores[0].end(), 0.0) / pvalores[0].size();			
				// pvalues2Val= std::accumulate( pvalores[1].begin(), pvalores[1].end(), 0.0) / pvalores[1].size();
				// pvalues3Val= std::accumulate( pvalores[2].begin(), pvalores[2].end(), 0.0) / pvalores[2].size();
				// pvalues4Val= std::accumulate( pvalores[3].begin(), pvalores[3].end(), 0.0) / pvalores[3].size();
				// pvaluesbis1Val= std::accumulate( pvaloresbis[0].begin(), pvaloresbis[0].end(), 0.0) / pvaloresbis[0].size();			
				// pvaluesbis2Val= std::accumulate( pvaloresbis[1].begin(), pvaloresbis[1].end(), 0.0) / pvaloresbis[1].size();
				// pvaluesbis3Val= std::accumulate( pvaloresbis[2].begin(), pvaloresbis[2].end(), 0.0) / pvaloresbis[2].size();
				// pvaluesbis4Val= std::accumulate( pvaloresbis[3].begin(), pvaloresbis[3].end(), 0.0) / pvaloresbis[3].size();
				// cout << pvalues1Val << endl;
				// cout << pvalues2Val << endl;
				// cout << pvalues3Val << endl;
				// cout << pvalues4Val << endl;
				//largest diff
				
				


				outputFile->cd();
				
				// pvaluesTree->Fill();

				pvalores[0].clear(); pvalores[1].clear(); pvalores[2].clear(); pvalores[3].clear();
				pvaloresbis[0].clear(); pvaloresbis[1].clear(); pvaloresbis[2].clear(); pvaloresbis[3].clear();


			
			}

			////////////////////////////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////////////////////////////

			
			count(largestdiff1Val, largest_diff_sim, 0);
			count(largestdiff2Val, largest_diff_sim, 1);
			count(largestdiff3Val, largest_diff_sim, 2);
			count(largestdiff4Val, largest_diff_sim, 3);

			largestdiffsTree->Fill();
			largest_diff_sim.clear();

			
		}

	}
	// pvaluesTree->Write();
	largestdiffsTree->Write();
	outputFile->Write();
	outputFile->Close();

	t.Stop();
	t.Print();

	return 0;
}


