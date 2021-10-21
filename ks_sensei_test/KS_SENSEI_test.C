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


//// Selection options ////

int x; 
int y; 

vector<vector<double>> pvalores(4);
vector<vector<double>> pvaloresbis(4);
vector<vector<double>> largest_diff_sim(4);


//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////


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
	// cout << largest_diff_sim[index1].size() << endl;
	// cout << "[";
	// cout << largest_diff_sim[index1][0] << ",";
	int counter=0;
	for (size_t i = 1; i < largest_diff_sim[index1].size(); i++)
	{
		if (largest_diff_sim[index1][0]>largest_diff_sim[index1][i])
		{
			largestdiffVal++;
		}
		if (largest_diff_sim[index1][0]==largest_diff_sim[index1][i])
		{
			// largestdiffVal=largestdiffVal+0.5;
			counter++;
		}
		
		// cout << largest_diff_sim[index1][i] <<",";
	}
	// cout << "]" << endl;
	largestdiffVal=largestdiffVal/(largest_diff_sim[index1].size()-counter);
	// cout << largestdiffVal << endl;
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


void analyse(vector<vector<double>> &events, int Entries, int hdu, int neighbours)
{

	auto test = TRandom3(0);
	test.SetSeed(0);

	for(int i_event=0;i_event<Entries; i_event++){
		double  xx= test.Uniform(0,1)*116000.0;
		events[0].push_back(xx);
	}
	test.SetSeed(0);
	for(int i_event=0;i_event<Entries; i_event++){
		double yy= test.Uniform(0,1)*116000.0;
		events[1].push_back(yy);
	}
	for(int i_event=0;i_event<Entries; i_event++){
			events[2].push_back(sqrt(pow(events[0][i_event],2)+pow(events[1][i_event],2)));
	}

	const int n = events[0].size();

	double dx, dy, distance;

	vector <float_t> distancetemp(neighbours,116000);
	for (int i = 0; i < n; ++i){
		distancetemp.assign(distancetemp.size(),116000);
		

		if (i>n-neighbours){distancetemp.pop_back();}
		// for (int j = i+1; j < n; ++j){
		for (int j = 0; j < n; ++j){
			if (i==j){continue;}
			dx=events[0][i]-events[0][j];
			dy=events[1][i]-events[1][j];
			distance=sqrt(pow(dx,2)+pow(dy,2));
			if (distance<distancetemp[distancetemp.size()-1])
			{
				distancetemp[distancetemp.size()-1]=distance;
				// cout << distance << endl;
				sort(distancetemp.begin(), distancetemp.end());
			}
		}

		for (size_t temp = 0; temp < distancetemp.size(); temp++)
		{
			// events[3].push_back(distancetemp[temp]);
			events[3].push_back(std::accumulate(distancetemp.begin(), distancetemp.end(), 0.0)/distancetemp.size());
		}
	}

	const int c = events[3].size();
}


void kolmogorovTest2(vector<vector<double>> ON, vector<vector<double>> OFF, vector<vector<double>> x_ON, vector<vector<double>> x_OFF, double &pvalue, int index1, int index2)
{

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

	double pvalue_x=0; double pvalue_y=0; double pvalue_c=0; double pvalue_d=0;

	// kolmogorovTest(ON, OFF, x_ON, x_OFF, pvalue_x,0,0);
	// kolmogorovTest(ON, OFF, x_ON, x_OFF, pvalue_y,1,0);
	// kolmogorovTest(ON, OFF, x_ON, x_OFF, pvalue_d,2,0);
	// kolmogorovTest(ON, OFF, x_ON, x_OFF, pvalue_c,3,1);
	kolmogorovTest2(ON, OFF, x_ON, x_OFF, pvalue_x,0,0);
	kolmogorovTest2(ON, OFF, x_ON, x_OFF, pvalue_y,1,0);
	kolmogorovTest2(ON, OFF, x_ON, x_OFF, pvalue_d,2,0);
	kolmogorovTest2(ON, OFF, x_ON, x_OFF, pvalue_c,3,1);

	return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[]){

	TStopwatch t; // time counter
	t.Start();

	gROOT->SetStyle("Plain");
	cout.precision(4);  // Print using four decimals precision
	cout<<fixed;        // ...exactly four decimals

	//outfilename
	std::string outfilename=("pvalues_"+std::to_string(atoi(argv[1]))+".root").c_str();

	//initialize output tree
	TFile *outputFile = new TFile(("files/"+outfilename).c_str(),"RECREATE"); //output file
	cout << outputFile->GetName() << endl;

	TTree * largestdiffsTree = new TTree("largestdiffs","largestdiffs");
	Float_t largestdiff1Val; Float_t largestdiff2Val; Float_t largestdiff3Val; Float_t largestdiff4Val;
	largestdiffsTree->Branch("largestdiff1",&largestdiff1Val,"largestdiff1/F");
	largestdiffsTree->Branch("largestdiff2",&largestdiff2Val,"largestdiff2/F");
	largestdiffsTree->Branch("largestdiff3",&largestdiff3Val,"largestdiff3/F");
	largestdiffsTree->Branch("largestdiff4",&largestdiff4Val,"largestdiff4/F");


	//set number of neighbours 

	// First set
	vector<vector<double>> THEO(4); //x, y , x**2+y**2, distance entre eventos
	int Entries = 10;
	int neighbours=Entries-1;
	analyse(THEO, Entries, 1, neighbours);

	//analyse and compare data
	for (int hdu = 1; hdu < 1000; hdu++){

		

		// Second set
		Entries = 10;
		vector<vector<double>> ON(4); //x, y , x**2+y**2, distance entre eventos
		analyse(ON, Entries, hdu , neighbours);

		// perform Kolmogorov test between the two sets
		kolmogorov(THEO, ON, hdu, outfilename);
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////

	count(largestdiff1Val, largest_diff_sim, 0);
	count(largestdiff2Val, largest_diff_sim, 1);
	count(largestdiff3Val, largest_diff_sim, 2);
	count(largestdiff4Val, largest_diff_sim, 3);
	outputFile->cd();
	largestdiffsTree->Fill();
	largestdiffsTree->Write();
	outputFile->Write();
	outputFile->Close();

	t.Stop();
	t.Print();

	return 0;
}


