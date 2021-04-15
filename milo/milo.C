// Autor: Mariano Cababie

#include <random>
#include <sys/time.h>
#include <time.h>
#include <cstdlib>
#include "fitsio.h"
#include <stdio.h>
#include <stdlib.h>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <iostream>
#include <chrono>
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TStopwatch.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include <string>
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TLeafC.h"
#include <sstream>
#include <algorithm>
#include <iterator>
#include <numeric>
#include "TRandom3.h"
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

//some global variables
// Char_t NROW;
int NROW;
int NCOL;
int NBINROW;
int NBINCOL;
int CCDNCOL;
int CCDNROW;
int CCDNPRES;
int RUNID;
int halo;
bool fake;


int compute1eEvents(TFile* &file, TFile* &outfile, Int_t N, vector<Int_t> &x1e, vector<Int_t> &y1e, vector<Double_t> &ePix1e, vector<Int_t> &ohdu1e, vector<Int_t> &mask1e, vector<Float_t> &distance1e, vector<Int_t> runID1e, const string &cutsPix, const string &cutsEmptyPix){


    TTree* dataPix = nullptr;	
	file->GetObject("calPixTree",dataPix);
    dataPix->SetName("dataPix");

    TTree* dataPixSel = nullptr; 
    dataPixSel = dataPix->CopyTree(cutsPix.c_str());
    dataPixSel->SetName("dataPixSel");
	dataPixSel->SetBranchStatus("*",1);

    TTreeReader readPix("dataPixSel", outfile);
	TTreeReaderValue<Int_t> x1electron(readPix, "x");
	TTreeReaderValue<Int_t> y1electron(readPix, "y");
	TTreeReaderValue<Double_t> epix1electron(readPix, "ePix");
	TTreeReaderValue<Int_t> ohdu(readPix, "ohdu");
	TTreeReaderValue<Int_t> mask(readPix, "mask");
    TTreeReaderValue<Float_t> distance1electron(readPix, "distance");
    
    TTree imgSumm("imgSumm","imgSumm");
    imgSumm.Branch("N",&N,"N/I");
    imgSumm.Branch("x1e",&x1e[0],"x1e[N]/I");
	imgSumm.Branch("y1e",&y1e[0],"y1e[N]/I");
	imgSumm.Branch("ePix1e",&ePix1e[0],"ePix1e[N]/D");
	imgSumm.Branch("ohdu1e",&ohdu1e[0],"ohdu1e[N]/I");
	imgSumm.Branch("mask1e",&mask1e[0],"mask1e[N]/I");
    imgSumm.Branch("runID1e",&runID1e[0],"runID1e[N]/I");
    imgSumm.Branch("distance1e",&distance1e[0],"distance1e[N]/F");


	while (readPix.Next()) {
		x1e.push_back(*x1electron);
		y1e.push_back(*y1electron);
		ohdu1e.push_back(*ohdu);
		mask1e.push_back(*mask);
        runID1e.push_back(RUNID);
		ePix1e.push_back(*epix1electron);
        distance1e.push_back(*distance1electron);
	}

    N=x1e.size();


    int a=0;
    vector<Int_t> xfake; vector<Int_t> yfake; vector<Int_t> ohdufake; vector<Int_t> maskfake;
    imgSumm.SetBranchAddress("x1e",&x1e[0]);
	imgSumm.SetBranchAddress("y1e",&y1e[0]);
    imgSumm.SetBranchAddress("ePix1e",&ePix1e[0]);
    imgSumm.SetBranchAddress("ohdu1e",&ohdu1e[0]);
    imgSumm.SetBranchAddress("mask1e",&mask1e[0]);
    imgSumm.SetBranchAddress("runID1e",&runID1e[0]);
    imgSumm.SetBranchAddress("distance1e",&distance1e[0]);
    imgSumm.Fill();
    imgSumm.Write();
    
    delete dataPix;
    delete dataPixSel;
    return 0;
}



int createOtherTree(TFile* &file, TFile* &outfile, const string &cutsPix, const string &cutsEmptyPix){

    TTree* dataPix = nullptr;	
	file->GetObject("calPixTree",dataPix);
    dataPix->SetName("dataPix");
    

    TTree* otherTreeAux = nullptr; 
    otherTreeAux = dataPix->CopyTree(cutsEmptyPix.c_str());
    otherTreeAux->SetName("otherTreeAux");
	otherTreeAux->SetBranchStatus("*",1);
    TTreeReader readOtherTreeAux("otherTreeAux", outfile);
    TTreeReaderValue<Int_t> maskAux(readOtherTreeAux, "mask");
    TTreeReaderValue<Int_t> ohduAux(readOtherTreeAux, "ohdu");
    TTreeReaderValue<Float_t> distanceAux(readOtherTreeAux, "distance");
    vector<Int_t>    masksOther;
    vector<Int_t>    ohdusOther;
    vector<Float_t>    distancesOther;
    vector<Int_t>    runIDOther;
    int N;

    TTree otherTree("otherTree","otherTree");
    otherTree.Branch("N",&N,"N/I");
    otherTree.Branch("mask",&masksOther[0],"mask[N]/I");
    otherTree.Branch("ohdu",&ohdusOther[0],"ohdu[N]/I");
    otherTree.Branch("distance",&distancesOther[0],"distance[N]/F");
    otherTree.Branch("runID",&runIDOther[0],"runID[N]/I");

    while (readOtherTreeAux.Next()) {
        masksOther.push_back(*maskAux);
        ohdusOther.push_back(*ohduAux);
		distancesOther.push_back(*distanceAux);
        runIDOther.push_back(RUNID);
    }

    N=masksOther.size();

    otherTree.SetBranchAddress("mask",&masksOther[0]);
    otherTree.SetBranchAddress("ohdu",&ohdusOther[0]);
    otherTree.SetBranchAddress("distance",&distancesOther[0]);
    otherTree.SetBranchAddress("runID",&runIDOther[0]);

    otherTree.Fill();
    otherTree.Write();

    delete otherTreeAux;
    delete dataPix;

    return 0; 
}


int analyseImage(const string &infile, const int &halo){

    int status = 0;

    TFile *file = TFile::Open(infile.c_str());
	TTree* data = nullptr;
	file->GetObject("hitSumm",data);

    //set file name
	std::vector<std::string> words;split(infile,words,'/');
    std::string outfilename;
    if (fake==true){
	    outfilename="fakemilo_"+words.back();
    }else{
        outfilename="milo_"+words.back();
    }
	
    TFile *outfile = new TFile(outfilename.c_str(),"RECREATE");

    cout << "\nProcessing file " << words.back() << endl;

	//getting NCOL and NROW
	TTree* config =nullptr; file->GetObject("headerTree_1",config); config->GetEvent(0);
	
	char* NROWchar;NROWchar=((TLeafC *) config->GetBranch("NROW")->GetLeaf("string"))->GetValueString();istringstream buffer(NROWchar); buffer >> NROW; 
	char* NCOLchar;NCOLchar=((TLeafC *) config->GetBranch("NCOL")->GetLeaf("string"))->GetValueString();istringstream buffer2(NCOLchar); buffer2 >> NCOL; 
	char* CCDNROWchar;CCDNROWchar=((TLeafC *) config->GetBranch("CCDNROW")->GetLeaf("string"))->GetValueString();istringstream buffer3(CCDNROWchar); buffer3 >> CCDNROW; 
	char* CCDNCOLchar;CCDNCOLchar=((TLeafC *) config->GetBranch("CCDNCOL")->GetLeaf("string"))->GetValueString();istringstream buffer4(CCDNCOLchar); buffer4 >> CCDNCOL; 
	char* CCDNPRESchar;CCDNPRESchar=((TLeafC *) config->GetBranch("CCDNPRES")->GetLeaf("string"))->GetValueString();istringstream buffer5(CCDNPRESchar); buffer5 >> CCDNPRES; 
    char* NBINROWchar;NBINROWchar=((TLeafC *) config->GetBranch("NBINROW")->GetLeaf("string"))->GetValueString();istringstream buffer6(NBINROWchar); buffer6 >> NBINROW; 
	char* NBINCOLchar;NBINCOLchar=((TLeafC *) config->GetBranch("NBINCOL")->GetLeaf("string"))->GetValueString();istringstream buffer7(NBINCOLchar); buffer7 >> NBINCOL; 
    // NBINCOL=1;
    // NBINROW=1;
    char* RUNIDchar;RUNIDchar=((TLeafC *) config->GetBranch("RUNID")->GetLeaf("string"))->GetValueString();istringstream buffer8(RUNIDchar); buffer8 >> RUNID; 


    
    std::string cutsPix="ePix>0.63 && ePix<2.5 && !(mask & (1+4+16+128+256+512+1024+4096)) && x>"+std::to_string(CCDNPRES+1)+" && x<="+std::to_string((CCDNCOL/2+CCDNPRES)/NBINCOL)+" && y>0 && y<"+std::to_string(min(NROW,CCDNROW/(2*NBINROW)))+" && (ohdu==4 || ohdu==3 || ohdu==2 || ohdu==1)";

    std::string cutsEmptyPix="!(mask & (1+4+16+128+256+512+1024+4096)) && x>"+std::to_string(CCDNPRES+1)+" && x<="+std::to_string((CCDNCOL/2+CCDNPRES)/NBINCOL)+" && y>0 && y<"+std::to_string(min(NROW,CCDNROW/(2*NBINROW)))+" && (ohdu==4 || ohdu==3 || ohdu==2 || ohdu==1)";
    
    // Clone Tree
	// TTree* hitSumm = nullptr;
	// hitSumm = data->CopyTree("");
	// hitSumm->SetName("hitSumm");
	// hitSumm->SetBranchStatus("*",1);

    // Clone configTree
	TTree* headerTree_1 = nullptr;
	headerTree_1 = config->CopyTree("");
	headerTree_1->SetName("headerTree_1");
	headerTree_1->SetBranchStatus("*",1);

	// Select 1e- events
	Int_t N;
    vector<Int_t>    x1e;
	vector<Int_t>    y1e;
	vector<Double_t>    ePix1e;
	vector<Int_t>    ohdu1e;
	vector<Int_t>    mask1e;
    vector<Float_t>    distance1e;
    vector<Int_t>    runID1e;


    cout << "Cut for 1e events = " << cutsPix <<"\n"<<endl;

    compute1eEvents(file, outfile, N, x1e, y1e, ePix1e, ohdu1e,mask1e, distance1e, runID1e, cutsPix, cutsEmptyPix);

    createOtherTree(file,outfile, cutsPix, cutsEmptyPix);

    
    // delete hitSumm;
    outfile->Write();
	outfile->Close();

	return status;
}

int main(int argc, char* argv[]){

time_t start,end;
double dif;
time (&start);

std::string infile(argv[1]);
halo=60; //global variable
fake=false;

int status = analyseImage(infile, halo);

time (&end);
dif = difftime (end,start);
cout << "\n\n---------------------------------" <<endl;
cout << "\n\nTook me " << dif << " seconds.\n\n" << endl;
cout << "---------------------------------\n\n" <<endl;
return 0;

}