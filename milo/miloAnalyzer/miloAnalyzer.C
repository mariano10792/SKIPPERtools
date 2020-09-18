#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include "TFile.h"
#include "TObject.h"
#include "TKey.h"
#include "TString.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TGraphAsymmErrors.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TLine.h"
#include "TLegend.h"
#include <TROOT.h>
#include "TSystem.h"
#include "TGaxis.h"
#include <cstdlib>
#include <string>
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TObjString.h"
#include "TH2D.h"
#include <fstream>
#include "TRandom3.h"
#include <iterator>




// Int_t N;
// vector<Int_t>    x1e;
// vector<Int_t>    y1e;
// vector<Double_t>    ePix1e;
// vector<Int_t>    ohdu1e;
// vector<Int_t>    mask1e;
// vector<Float_t>    distance1e;
// Int_t    runID1e;

using namespace std;

// void Enable_and_Set_Branches(TTree* & tree);



void miloAnalyzer(){



  // TFile *file = TFile::Open("/mnt/mariano/MINOSnewfirmware/run3_NROW500_NBINROW2_NCOL240_NBINCOL2_EXPOSURE0/0/poisson.root");

  TFile *file = TFile::Open("/mnt/mariano/SNOLAB/0/poisson.root");
  
  TTree *tree = (TTree*) file->Get("imgSumm");
  TTree *otherTree = (TTree*) file->Get("otherTree");
  // Enable_and_Set_Branches(tree);
  tree->SetBranchStatus("*",1); //enable all branches
  otherTree->SetBranchStatus("*",1); //enable all branches
	

  TTreeReader readPix("imgSumm", file);
	TTreeReaderArray<Int_t> x1e(readPix, "x1e");
	TTreeReaderArray<Int_t> y1e(readPix, "y1e");
	TTreeReaderArray<Double_t> ePix1e(readPix, "ePix1e");
	TTreeReaderArray<Int_t> ohdu1e(readPix, "ohdu1e");
	TTreeReaderArray<Int_t> mask1e(readPix, "mask1e");
  TTreeReaderArray<Float_t> distance1e(readPix, "distance1e");

  TTreeReader readOtherTree("otherTree", file);
  TTreeReaderArray<Int_t> otherMask(readOtherTree, "mask");
  TTreeReaderArray<Float_t> otherDistance(readOtherTree, "distance");

  // auxilar variables
  vector<Float_t> Nvector;
  vector<Float_t> rate;
  vector<Float_t> rateSim;
  int j=0; // some counter

  // setting plots
  int bins =50;



  //cuts
  int distance=100;
  int xmin=8;
  int xmax=226;
  int ymin=0;
  int ymax=501;

  

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////



  TH1D * h  =  new TH1D("1eEvents","",bins,0,50);
  TCanvas * c1   = new TCanvas("Spectrum","Spectrum",1920,1080);
  c1->Divide(2,3);
  c1->cd(1);

  Nvector.clear();
  while (readPix.Next()) {
    int N=0;
    for (int i = 0, ni =  x1e.GetSize(); i < ni; ++i) {
        if (distance1e[i]>distance && x1e[i]>xmin && x1e[i]<xmax && y1e[i]>ymin && y1e[i]<ymax && ePix1e[i]<3 && (ohdu1e[i]==1 || ohdu1e[i]==2)){
          N++;
        }
    }
    Nvector.push_back(N);

    h->Fill(N);
    h->SetTitle(("distance<"+std::to_string(distance)+" && xmin>"+std::to_string(xmin)+" && xmax<"+std::to_string(xmax)+" && ymin>"+std::to_string(ymin)+"&& ymax<"+std::to_string(ymax)).c_str());
  }
  h->Draw();

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////


  TH1D * h2  =  new TH1D("1eEvents(unnormalized)","",bins,0,.0002);
  c1->cd(2);
  rate.clear();
  while (readOtherTree.Next()) {
    int emptyN=0;
    for (int i = 0, ni =  otherDistance.GetSize(); i < ni; ++i) {
      if (otherDistance[i]>distance){
        emptyN++;
      }     
    }
    rate.push_back(Nvector[j]/emptyN);
    h2->Fill(Nvector[j]/emptyN);
    j++;
  }
  h2->Draw();
  
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////




  
  TH1D * h3  =  new TH1D("1eEvents(normalized)","",bins,0,45);
  c1->cd(3);
  readOtherTree.Restart();
  j=0;
  vector <Float_t> counts;
  while (readOtherTree.Next()) {
    int emptyN=0;
    for (int i = 0, ni =  otherDistance.GetSize(); i < ni; ++i) {
      if (otherDistance[i]>distance){
        emptyN++;
      }
      
    }
    h3->Fill((Nvector[j]*221500)/emptyN);
    counts.push_back((Nvector[j]*221500)/emptyN);
    j++;
  }
  h3->Draw();

  //////////////////////////////////////// Kolmogorov //////////////////////////////////////// 

  std::sort(counts.begin(),counts.end());
  for(float i : counts)
    cout << i << endl;
    std::ofstream output_file("./counts.txt");
  copy(counts.begin(), counts.end(), ostream_iterator<float>(output_file , ",")); 



  //////////////////////////////////////// ChiSquare //////////////////////////////////////// 
  // TRandom3 *random=new TRandom3(0);
  
  // Float_t mean = std::accumulate(rate.begin(),rate.end(),0.00000000)/rate.size();
  // vector <Float_t> chi2;  
  // vector <Float_t> chi2Sim;  
  // for(float i : rate){
  //   chi2.push_back(pow(i-mean,2)/mean);
  //   // cout << random->Poisson(mean*h3->GetMean())/h3->GetMean() << endl;
  //   // cout << i << endl;
  //   chi2Sim.push_back(pow(random->Poisson(mean*h3->GetMean())/h3->GetMean()-mean,2)/mean);
  // }
  // cout << (std::accumulate(chi2.begin(),chi2.end(),0.00000000)/chi2.size())*h3->GetMean() << endl;
  // cout << (std::accumulate(chi2Sim.begin(),chi2Sim.end(),0.00000000)/chi2Sim.size())*h3->GetMean() << endl;



  //////////////////////////////////////// ChiSquare (gaussian)//////////////////////////////////////// 

  Float_t mean = std::accumulate(counts.begin(),counts.end(),0.00000000)/counts.size();
  vector <Float_t> chi2;
  for(float i : counts){
    chi2.push_back(pow(i-mean,2)/mean);
  }

  cout << "chi2 (gaus) = "  <<std::accumulate(chi2.begin(),chi2.end(),0.00000000)/(chi2.size()-1) << endl;
  

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////


  //new cuts
  distance=60;
  TH1D * h4  =  new TH1D("Electrons4","",bins,0,50);
  c1->cd(4);
  
  readPix.Restart();
  Nvector.clear();
  while (readPix.Next()) {
    int N=0;
    for (int i = 0, ni =  x1e.GetSize(); i < ni; ++i) {
        if (distance1e[i]>distance && x1e[i]>xmin && x1e[i]<xmax && y1e[i]>ymin && y1e[i]<ymax && ePix1e[i]<3 && (ohdu1e[i]==1 || ohdu1e[i]==2)){
          N++;
        }
    }
    Nvector.push_back(N);
    h4->Fill(N);
    h4->SetTitle(("distance<"+std::to_string(distance)+" && xmin>"+std::to_string(xmin)+" && xmax<"+std::to_string(xmax)+" && ymin>"+std::to_string(ymin)+"&& ymax<"+std::to_string(ymax)).c_str());
  }
  h4->Draw();

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////

  TH1D * h5  =  new TH1D("Electrons5","",bins,0,.0002);
  c1->cd(5);
  rate.clear();
  j=0;
  readOtherTree.Restart();
  while (readOtherTree.Next()) {
    int emptyN=0;
    for (int i = 0, ni =  otherDistance.GetSize(); i < ni; ++i) {
      if (otherDistance[i]>distance){
        emptyN++;
      }
      
    }
    rate.push_back(Nvector[j]/emptyN);
    h5->Fill(Nvector[j]/emptyN);
    j++;
  }
  h5->Draw();

  TH1D * h6  =  new TH1D("Electrons6","",bins,0,250000);
  c1->cd(6);
  readOtherTree.Restart();
  while (readOtherTree.Next()) {
    int emptyN=0;
    for (int i = 0, ni =  otherDistance.GetSize(); i < ni; ++i) {
      if (otherDistance[i]>distance){
        emptyN++;
      }
      
    }
    h6->Fill(emptyN);
    j++;
  }
  h6->Draw();

  mean = std::accumulate(rate.begin(),rate.end(),0.00000000)/rate.size();
  chi2.clear();  
  for(float i : rate)
    chi2.push_back(pow(i-mean,2)/mean);

  cout << (std::accumulate(chi2.begin(),chi2.end(),0.00000000))*h6->GetMean() << endl;
  cout << (std::accumulate(chi2.begin(),chi2.end(),0.00000000)/chi2.size())*h6->GetMean() << endl;



  // c1->cd(4);
  // tree->Draw("x1e:y1e",("distance1e<"+std::to_string(distance)+" && x1e>"+std::to_string(xmin)+" && x1e<"+std::to_string(xmax)+" && y1e>"+std::to_string(ymin)+"&& y1e<"+std::to_string(ymax)).c_str(),"colz");


}




