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



using namespace std;

int kmin=5;
int kmax=90;
// int kmin=50;
// int kmax=120;

void otherLikelihoodTest(vector<Float_t> &sim, Float_t &Ratio){
  vector <Float_t> meanaux;
  vector <Float_t> ratio;
  for (int i = 0; i < sim.size(); ++i){
    meanaux.push_back(sim[i]);
  }
  for (int i = 0; i < sim.size(); ++i){
    ratio.push_back(2*sim[i]*log(sim[i]/(std::accumulate(meanaux.begin(), meanaux.end(), 0.00000)/meanaux.size())));
  }
  Ratio = std::accumulate(ratio.begin(), ratio.end(), 0.00000);
}

void CramerVonMisesTests(vector<Float_t> &sim, TH1D* &h3, vector <Float_t> &CkAux , vector <Float_t> &A2VECTOR){
  
  ///// First test /////
  vector <Float_t> Aj;
  vector <Float_t> Bj;
  vector <Float_t> Ck;    
  Float_t AjAux;
  Float_t BjAux;
  Float_t CjAux;

  ///// Second test /////
  vector <Float_t> obs;
  vector <Float_t> pj;
  vector <Float_t> Hj;
  vector <Float_t> zj;
  vector <Float_t> Zj;
  vector <Float_t> A2aux;
  for (int k = kmin, N = kmax; k < N; ++k){
    int K=k-1; //bins start at 1 so the bin of n==0 is the k==1 bin.
    
    ///// First test /////

    Aj.push_back(h3->GetBinContent(k));
    
    Bj.push_back(exp(-h3->GetMean())*pow(h3->GetMean(),K)/tgamma(K+1));
    
    AjAux = (1/h3->GetEntries())*std::accumulate(Aj.begin(), Aj.end(),0.00000);
    
    BjAux = std::accumulate(Bj.begin(), Bj.end(),0.00000);
    
    CjAux = (1/h3->GetEntries())*h3->GetBinContent(k);
    
    Ck.push_back(h3->GetEntries()*pow(AjAux-BjAux,2) * CjAux);
    
    // cout << h3->GetEntries()*pow(AjAux-BjAux,2) * CjAux << " | " << h3->GetBinContent(k) << " | " << K <<endl;
    // cout <<   AjAux << " | " << BjAux <<  " | " << h3->GetBinContent(k) << " | " << K <<endl;
    
    ///// Second test ////
    obs.push_back(h3->GetBinContent(k));
    pj.push_back(exp(-h3->GetMean())*pow(h3->GetMean(),K)/tgamma(K+1));
    zj.push_back(obs[K]-(h3->GetEntries()*pj[K]));
    Zj.push_back(std::accumulate(zj.begin(), zj.end(), 0.00000));
    Hj.push_back(std::accumulate(pj.begin(), pj.end(), 0.00000));
    A2aux.push_back((pow(Zj[K],2)*pj[K])/(h3->GetEntries()*Hj[K]*(1-Hj[K])));
    
    

  }
  A2VECTOR.push_back(std::accumulate(A2aux.begin(), A2aux.end(), 0.00000));
  CkAux.push_back(std::accumulate(Ck.begin(), Ck.end(), 0.00000));
  
}



void miloAnalyzer(){


  string fileName = "/mnt/mariano/MINOSnewfirmware/run3_NROW500_NBINROW2_NCOL240_NBINCOL2_EXPOSURE0/0/milo.root";
  // string fileName = "/mnt/mariano/SNOLAB/0/milo.root";
  // string fileName = "/mnt/mariano/beforeSENSEI2020/shodata/halo60/milo.root";
  
  TFile *file = TFile::Open(fileName.c_str());
  // TFile *file = TFile::Open("/mnt/mariano/SNOLAB/0/milo.root");

  // setting plots
  int binmin =0;
  int binmax =40;
  int bins =(binmax-binmin);



  //cuts
  vector<Int_t> ohdus={1,1,2,2};
  // vector<Int_t> ohdus={2,2,3,3};
  // int distance=60;
  int distance=100;
  int xmin=8;
  int xmax=226;
  // int xmax=450;
  int ymin=0;
  int ymax=501;
  // int ymax=3072;


  
  cout << "Processing " << fileName << endl;
  std::string outfilename="output.root";
	TFile *outfile = new TFile(outfilename.c_str(),"RECREATE");
  
  TTree *tree = (TTree*) file->Get("imgSumm");
  TTree *otherTree = (TTree*) file->Get("otherTree");
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
  vector<Int_t> emptyNvector;
  vector<Double_t> Nvector; //need double for likelihood!
  vector<Double_t> Nvectorobs; //need double for likelihood!
  vector<Float_t> rate;
  vector<Float_t> rateSim;
  int j=0; // some counter
  vector<Int_t> excluded;
  int Nmin=0; int Nmax=19;
  int edgecut=4; //64






  ///////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////

  TCanvas * c1   = new TCanvas("Spectrum","Spectrum",1920,1080);
  c1->Divide(2,3);


  TH1D * h  =  new TH1D("Events1e","",bins,binmin,binmax);
  c1->cd(1);
  Nvector.clear();
  j=0;
  while (readPix.Next()) {
    int N=0;
    for (int i = 0, ni =  x1e.GetSize(); i < ni; ++i) {
        if (distance1e[i]>distance && !(mask1e[i] & (edgecut)) &&  x1e[i]>xmin && x1e[i]<xmax && y1e[i]>ymin && y1e[i]<ymax && ePix1e[i]<3 && (ohdu1e[i]==ohdus[0] || ohdu1e[i]==ohdus[1] || ohdu1e[i]==ohdus[2] || ohdu1e[i]==ohdus[3])){
          N++;
        }
    }
    Nvector.push_back(N);
    if (N>Nmax || N<Nmin){
      excluded.push_back(j);
    }else{     
      // cout << N << endl;
      Nvectorobs.push_back(N);
      h->Fill(N);
    }
    j++;
  }
  h->SetTitle(("distance>"+std::to_string(distance)+" && xmin>"+std::to_string(xmin)+" && xmax<"+std::to_string(xmax)+" && ymin>"+std::to_string(ymin)+"&& ymax<"+std::to_string(ymax)).c_str());
  h->Draw();






  ///////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////


  TH1D * h2  =  new TH1D("Rates1e.","",bins,0,.0006);
  c1->cd(2);
  rate.clear();
  j=0;
  bool go; 
  while (readOtherTree.Next()) {
    int emptyN=0; go=true;
    for(float i : excluded){if(i==j){cout << "Excluded "<< Nvector[j] << endl; go=false;}}
    if (go==false){j++;continue;}
    for (int i = 0, ni =  otherDistance.GetSize(); i < ni; ++i) {
      if (otherDistance[i]>distance && !(otherMask[i] & (edgecut))){
        emptyN++;
      }     
    }
    emptyNvector.push_back(emptyN);
    rate.push_back(Nvector[j]/emptyN);
    h2->Fill(Nvector[j]/emptyN);
    j++;
  }
  h2->Draw();
  




  ///////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////

  
  TH1D * h3  =  new TH1D("Events1eNorm","",bins,binmin,binmax);
  c1->cd(3);
  readOtherTree.Restart();
  j=0;
  vector <Float_t> counts;
  vector <Int_t> countsInt;
  go=true;
  while (readOtherTree.Next()) {
    int emptyN=0; go=true;
    for(float i : excluded){if(i==j){cout << "Excluded "<< Nvector[j] << endl; go=false;}}
    if (go==false){j++;continue;}
    for (int i = 0, ni =  otherDistance.GetSize(); i < ni; ++i) {
      if (otherDistance[i]>distance && !(otherMask[i] & (edgecut))){
        emptyN++;
      }     
    }  
    counts.push_back((Nvector[j]*2*(ymax-ymin)*(xmax-xmin))/emptyN);
    countsInt.push_back(round(counts[j]));
    h3->Fill((Nvector[j]*2*(ymax-ymin)*(xmax-xmin))/emptyN); //floats
    // h3->Fill(round(counts[j])); //rounded integers
    j++;
  }
  h3->Draw("HIST");
  
  

  



  
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  

  //////////////////////////////////////// Likelihood test ////////////////////////////////////////


  //initialize
  gRandom->SetSeed(0);
  Float_t mu=h2->GetMean(); // mu
  cout << "mu =" << mu << endl;  
  Int_t NSIMS=100000;
  vector<vector<Int_t>> mat;
  mat.resize(NSIMS, std::vector<Int_t>(h3->GetEntries()));

  // simulate events
  for (size_t i = 0; i < h3->GetEntries(); i++)
  {
    TF1 fPoissonBin("f",("ROOT::Math::poisson_pdf(x,"+std::to_string(mu*emptyNvector[i])+")").c_str(),0,200);
    for (size_t j = 0; j < NSIMS; j++)
    {
      mat[j][i]=fPoissonBin.GetRandom(); // each column is the same bin
    }
    
  }
  
  // calculate likelihoods
  vector<Float_t> Lj;
  vector<Float_t> pi;
	for (vector<Int_t> row: mat) {
    pi.clear();
    int i=0;
    for (Int_t value : row){
      pi.push_back(-log(exp(-mu*emptyNvector[i])*pow(mu*emptyNvector[i],value)/tgamma(value+1)));
      // cout << value <<endl;
      i++;
    }
    // cout << " | "<<endl;
    Lj.push_back(std::accumulate(pi.begin(), pi.end(), 0.00000));
    // cout << std::accumulate(pi.begin(), pi.end(), 0.00000) << endl;
	}
  cout << "Mean value: "<< std::accumulate(Lj.begin(), Lj.end(), 0.00000)/Lj.size() << endl;
  double max = *max_element(Lj.begin(), Lj.end());
  cout<<"Max value : "<<max<<endl;
  
  // calculate measured likelihood
  pi.clear();
  for (size_t i = 0; i < emptyNvector.size(); i++)
  {
    
    // pi.push_back(-log(exp(-rate[i]*emptyNvector[i])*pow(rate[i]*emptyNvector[i],Nvectorobs[i])/tgamma(Nvectorobs[i]+1)));
    pi.push_back(-log(exp(-mu*emptyNvector[i])*pow(mu*emptyNvector[i],Nvectorobs[i])/tgamma(Nvectorobs[i]+1)));
    
    // cout << rate[i]*emptyNvector[i] << " | " <<" | " << pi[i] << endl<< endl;
    // cout << pi[i] << " | " << rate[i]*empt yNvector[i] << " | "<< Nvectorobs[i] << " | " << i <<  endl<< endl;
  }
  

  Double_t pi_value= std::accumulate(pi.begin(), pi.end(), 0.00000);
  Float_t nLikelihood=0;
  for(float i : Lj) if(i>=pi_value){nLikelihood++;}
  cout << "p-value Likelihood = " << nLikelihood/NSIMS << endl;
  cout << nLikelihood << " | "  << NSIMS << endl; 
  cout << std::accumulate(pi.begin(), pi.end(), 0.00000) << endl;










  
  //////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////// Other tests //////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////
  
  
  TF1 fPoisson("f",("ROOT::Math::poisson_pdf(x,"+std::to_string(h3->GetMean())+")").c_str(),0,200);
  vector<Float_t> RATIO; //likelihood test
  vector <Float_t> CkAux; // Cramer vos Mises 1st test
  vector<Float_t> A2VECTOR; // Cramer vos Mises 2nd test
  vector<Float_t> sim;
  Int_t Nsims=1000;
  
  for (size_t j = 0; j < Nsims; j++)
  {
    

    ////////////////////////////// Simulations //////////////////////////////////////// 
    
    sim.clear();
    if(j==Nsims-1){
      for (int i = 0; i < counts.size(); ++i){
      sim.push_back(counts[i]);
      // sim.push_back(fPoisson.GetRandom());
      }
    }else{
      for (int i = 0; i < counts.size(); ++i){
      sim.push_back(fPoisson.GetRandom());
      }
    }



    //reset histogram because we are going to fill it with simulations
    h3->Reset("ICESM");
    for (int i = 0; i < sim.size(); ++i){h3->Fill(sim[i]);}
    h3->Draw();
       
       

    
    
    
    //////////////////////////////////////// Other test (likelihood)////////////////////////////////////////   
    Float_t Ratio=0;
    // otherLikelihoodTest(sim,Ratio);
    RATIO.push_back(Ratio);


    ///////////////////////////////////////////////////////////////////////////////////
    //////////// Cramer vos Mises test from  10.1016/S0378-3758(00)00114-2 ////////////
    ///////////////////////// Cramer vos Mises2 10.2307/3315735 ///////////////////////
    ///////////////////////////////////////////////////////////////////////////////////
    
    CramerVonMisesTests(sim, h3, CkAux, A2VECTOR);

    
    

  }
  
  // cout << "A2 mean = " << std::accumulate(A2VECTOR.begin(), A2VECTOR.end(), 0.00000)/A2VECTOR.size() << endl;
  cout << "CVM mean = "<< std::accumulate(CkAux.begin(), CkAux.end(), 0.00000)/CkAux.size() << endl;
  Float_t CVMvalue=CkAux[Nsims-1];
  cout << "CVM value = "<< CVMvalue << endl;
  std::sort (CkAux.begin(), CkAux.end());
  
  Float_t prob=-1; // by construction, it counts the observed event
  for(float i : CkAux) if(i>=CVMvalue){prob++;}
  cout << prob << " | "  << Nsims << endl;
  cout << "p-value CVM = " <<prob/(CkAux.size()-1) << endl;


  // h3->Rebin();

  ///////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////// Fit ///////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////

  
  // TF1 *f1 = new TF1("f1","[0]*TMath::Power(([1]/[2]),(x/[2]))*(TMath::Exp(-([1]/[2])))/TMath::Gamma((x/[2])+1.)", kmin, kmax); // "xmin" = 0, "xmax" = 10
  // f1->SetParameters(h3->GetEntries(), h3->GetMean(), 1); // you MUST set non-zero initial values for parameters
  // h3->Fit("f1", "R"); // "R" = fit between "xmin" and "xmax" of the "f1"
  // TF1 *fit = h3->GetFunction("f1");
  // cout << "ChiSquare = "<< fit->GetChisquare() << endl;
  // cout << "NDF = "<< fit->GetNDF() << endl;
  // cout << "P-value (Chi2) = "<< fit->GetProb() << endl;






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

  // Float_t meanDist = std::accumulate(counts.begin(),counts.end(),0.00000000)/counts.size();
  // vector <Float_t> chi2;
  // for(float i : counts){
  //   chi2.push_back(pow(i-meanDist,2)/meanDist);
  // }

  // cout << "chi2 (gaus) = "  <<std::accumulate(chi2.begin(),chi2.end(),0.00000000)/(chi2.size()-1) << endl;
  

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////

  

  // TH1D * h4  =  new TH1D("Unnamed","",bins,0,60);
  c1->cd(4);
  tree->Draw(("x1e:y1e>>h4("+std::to_string(bins/10)+",0,500,"+std::to_string(bins/10)+",0,230)").c_str(),("(ohdu1e=="+std::to_string(ohdus[0])+" || ohdu1e=="+std::to_string(ohdus[1])+" || ohdu1e=="+std::to_string(ohdus[2])+" || ohdu1e=="+std::to_string(ohdus[3])+") && distance1e>"+std::to_string(distance)+" && x1e>"+std::to_string(xmin)+" && x1e<"+std::to_string(xmax)+" && y1e>"+std::to_string(ymin)+"&& y1e<"+std::to_string(ymax)).c_str(),"colz");

  // readPix.Restart();
  // Nvector.clear();
  // while (readPix.Next()) {
  //   int N=0;
  //   for (int i = 0, ni =  x1e.GetSize(); i < ni; ++i) {
  //       if (distance1e[i]>distance && x1e[i]>xmin && x1e[i]<xmax && y1e[i]>ymin && y1e[i]<ymax && ePix1e[i]<3 && (ohdu1e[i]==1 || ohdu1e[i]==2)){
  //         N++;
  //       }
  //   }
  //   Nvector.push_back(N);
  //   h4->Fill(N);
  //   h4->SetTitle(("distance<"+std::to_string(distance)+" && xmin>"+std::to_string(xmin)+" && xmax<"+std::to_string(xmax)+" && ymin>"+std::to_string(ymin)+"&& ymax<"+std::to_string(ymax)).c_str());
  // }
  // h4->Draw();

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////

  c1->cd(5);
  
  TH1D * h5  =  new TH1D("Unnamed","",bins,0,.0002);
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



  ///////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////




  TH1D * h6  =  new TH1D("EMPTY_PIXELS","",bins,50000,250000);
  c1->cd(6);


  for (size_t i = 0; i < emptyNvector.size(); i++)
  {
    h6->Fill(emptyNvector[i]);
    // cout << emptyNvector[  i] << endl;
  }
  h6->Draw();


  // readOtherTree.Restart();
  // while (readOtherTree.Next()) {
  //   int emptyN=0;
  //   for (int i = 0, ni =  otherDistance.GetSize(); i < ni; ++i) {
  //     if (otherDistance[i]>distance){
  //       emptyN++;
  //     }
      
  //   }
  //   h6->Fill(emptyN);
  //   j++;
  // }
  // h6->Draw();

  // meanDist = std::accumulate(rate.begin(),rate.end(),0.00000000)/rate.size();
  // chi2.clear();  
  // for(float i : rate)
  //   chi2.push_back(pow(i-meanDist,2)/meanDist);

  // cout << (std::accumulate(chi2.begin(),chi2.end(),0.00000000))*h6->GetMean() << endl;
  // cout << (std::accumulate(chi2.begin(),chi2.end(),0.00000000)/chi2.size())*h6->GetMean() << endl;


  // outfile->Write();
	// outfile->Close();
  // c1->cd(3);
  // TF1 * fPoissonBin= new TF1("f",("1.8*ROOT::Math::poisson_pdf(x,"+std::to_string(h3->GetMean())+")").c_str(),0,200);
  // fPoissonBin->Draw("same");
  // h3->Scale(1./h3->GetEntries());
  // c1->Update();
}




