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
int halo;
bool dense;
bool old;

int processCommandLineArgs(const int argc, char *argv[], bool &dense, int &halo){

    int opt=0;
    while ( (opt = getopt(argc, argv, "h:do")) != -1) {
        switch (opt) {
            case 'd':
				dense=true;
				cout << "\n\nEnzo will use dense configuration. If you are seeing this then your file must have lots of High Energy Events." << endl;
                break;
			case 'h':
				cout << "\n\nArgument halo = " <<atoi(optarg) << endl;
				halo = atoi(optarg);
                break;
			case 'o': //old files don't have a binning branch so they crash.
				old=true;
				cout << "\n\nBinning set at 1." << endl;
                break;
            default: /* '?' */
                return 0;
        }
    }
    return 0;
}

int cutEvents(TFile* &file,TFile* &outfile, TTreeReader &readEvents, TTree* &hitSumm, vector<Int_t> &masks){

	readEvents.Restart();
	TTreeReaderValue<Int_t> xMin(readEvents, "xMin");
	TTreeReaderValue<Int_t> xMax(readEvents, "xMax");
	TTreeReaderArray<Int_t> xPix(readEvents, "xPix");
	TTreeReaderValue<Int_t> yMin(readEvents, "yMin");
	TTreeReaderValue<Int_t> yMax(readEvents, "yMax");
	TTreeReaderArray<Int_t> yPix(readEvents, "yPix");
	TTreeReaderValue<Int_t> ohdu(readEvents, "ohdu");
	TTreeReaderArray<Float_t> ePix(readEvents, "ePix");
	
	

	// I need a zombie tree so the Aux and the original don't point to the same address when their TTreeReaders are iterated (see Documentation of CopyTree)
	TTree* hitSummZombie = nullptr;
	hitSummZombie = hitSumm->CopyTree("");
	hitSummZombie->SetName("hitSummZombie");
	hitSummZombie->SetBranchStatus("*",1);

	TTree* hitSummAux = nullptr;
	hitSummAux = hitSummZombie->CopyTree("");
	hitSummAux->SetName("hitSummAux");
	hitSummAux->SetBranchStatus("*",1);

	TTreeReader compareEvents("hitSummAux", outfile);
	TTreeReaderArray<Int_t> xPix2(compareEvents, "xPix");
	TTreeReaderArray<Int_t> yPix2(compareEvents, "yPix");
	TTreeReaderArray<Float_t> ePix2(compareEvents, "ePix");
	TTreeReaderValue<Int_t> ohdu2(compareEvents, "ohdu");
	delete hitSummZombie;

	

	Int_t booleanVal;
	Int_t sum;
	TBranch *booleanBranch = hitSumm->Branch("boolean",&booleanVal,"boolean/I");
	TBranch *sumBranch = hitSumm->Branch("NEvent",&sum,"NEvent/I");
	int j1=0;
	


	std::string cutsEvents ="e>=100 && xVar>0 && yVar>0 && xMin>"+std::to_string(CCDNPRES)+" && xMax<="+std::to_string((CCDNCOL/2+CCDNPRES)/NBINCOL)+" && yMin>0 && yMax<"+std::to_string(min(NROW,CCDNROW/(2*NBINROW)))+" && ohdu<3";
	std::string cutsPix="ePix>0.63 && ePix<1.5 && !(mask & (1+4+16+128+512+1024)) && x>"+std::to_string(CCDNPRES)+" && x<="+std::to_string((CCDNCOL/2+CCDNPRES)/NBINCOL)+" && y>0 && y<"+std::to_string(min(NROW,CCDNROW/(2*NBINROW)))+" && ohdu<3";






	int counter=0;
	while (readEvents.Next()) {
		
		std::vector<int> area(NROW*NCOL*2, 0);
		compareEvents.Restart();
		hitSumm->GetEntry(j1);
		booleanVal=1;
		if (*xMin<CCDNPRES+halo || *xMax>CCDNCOL/2+CCDNPRES-halo || *yMin<halo || *yMax >min(NROW,CCDNROW/NBINROW)-halo){booleanVal=0; sum=0; booleanBranch->Fill(); sumBranch->Fill(); j1++; continue;}
		// if (*xMin<10+halo || *xMax>450-halo || *yMin<halo || *yMax >NROW-halo || !(*xMin>400+halo || *xMax<250-halo || *yMax<1000-halo)){booleanVal=0; booleanBranch->Fill(); j1++; continue;} #module4 on excesses
		
		for (int i = 0, ni =  xPix.GetSize(); i < ni; ++i) {
			if (booleanVal==0){cout << "I broke." << endl; break;}
			
			//halo radius for a given pixelmak
			int xStart = max(0,xPix[i]-halo);
			int xEnd = min(NCOL-1,xPix[i]+halo); 
			int yStart = max(0,yPix[i]-halo);
			int yEnd = min(NROW-1,yPix[i]+halo);


		
			for (int x = xStart; x <= xEnd; ++x){			
				for (int y = yStart; y <= yEnd; ++y){			
					if(!(masks[(x)+(y)*NCOL+NCOL*NROW*((*ohdu)-1)] & (1+2+4+16+128+512+1024))){
						area[(x)+(y)*NCOL+NCOL*NROW*((*ohdu)-1)]=1;

					}
				}
			}
		}
		sum=std::accumulate(area.begin(), area.end(), 0);		
		
		compareEvents.Restart();
		if (!dense){
			for (int i = 0, ni =  xPix.GetSize(); i < ni; ++i) {
				int j2=0;
				while (compareEvents.Next()){
					if (j1==j2 or *ohdu!=*ohdu2){j2++; continue;}
					for (int k = 0, nk =  xPix2.GetSize(); k < nk; ++k) {
						if (pow(pow(xPix[i]-xPix2[k],2)+pow(yPix[i]-yPix2[k],2),0.5)<2*halo && j1!=j2){
							booleanVal=0;
							break;
						}
					}		
					j2++;
					if (booleanVal==0){break;} //breaks the while
				}
				if (booleanVal==0){break;}
			}
		}






		booleanBranch->Fill();
		sumBranch->Fill();
		j1++;
	}
	delete hitSummAux;
	return 0;
}

int Compute1eEvents(TFile* &file, TFile* &outfile, vector<Int_t> &x1e, vector<Int_t> &y1e, vector<Double_t> &ePix1e, vector<Int_t> &ohdu1e, vector<Int_t> &mask1e, const string &cutsPix, vector<Int_t> &masks){

	TTree* dataPix = nullptr;	
	file->GetObject("calPixTree",dataPix);
	dataPix->SetName("dataPix");
	dataPix->SetBranchStatus("*",1);
	TTreeReader readCalPixTree("dataPix", file);
	TTreeReaderValue<Int_t> maskVal(readCalPixTree, "mask");
	TTreeReaderValue<Int_t> xVal(readCalPixTree, "x");
	TTreeReaderValue<Int_t> yVal(readCalPixTree, "y");
	TTreeReaderValue<Int_t> ohduVal(readCalPixTree, "ohdu");
	while (readCalPixTree.Next()) {
		if (*ohduVal<3){
			masks[(*xVal)+(*yVal)*NCOL+NCOL*NROW*((*ohduVal)-1)]=*maskVal;
		}
	}

	TTree* dataPixSel = nullptr;
	dataPixSel = dataPix->CopyTree(cutsPix.c_str());
	dataPixSel->SetName("dataPixSel");
	// cout << cutsPix << endl;
	dataPixSel->SetBranchStatus("*",1);

	TTreeReader readPix("dataPixSel", outfile);
	TTreeReaderValue<Int_t> x1electron(readPix, "x");
	TTreeReaderValue<Int_t> y1electron(readPix, "y");
	TTreeReaderValue<Double_t> epix1electron(readPix, "ePix");
	TTreeReaderValue<Int_t> ohdu(readPix, "ohdu");
	TTreeReaderValue<Int_t> mask(readPix, "mask");

	while (readPix.Next()) {
		x1e.push_back(*x1electron);
		y1e.push_back(*y1electron);
		ohdu1e.push_back(*ohdu);
		mask1e.push_back(*mask);
		ePix1e.push_back(*epix1electron);
	}

	delete dataPix;
	delete dataPixSel;
	return 0;	
}


int WriteTree(TTreeReader &readEvents,TTree* &hitSumm, vector<Int_t> &x1e, vector<Int_t> &y1e, vector<Double_t> &ePix1e, vector<Int_t> &ohdu1e, vector<Int_t> &mask1e){

	int j1=0;
	readEvents.Restart();
	TTreeReaderArray<Int_t> xPix(readEvents, "xPix");
	TTreeReaderArray<Int_t> yPix(readEvents, "yPix");
	TTreeReaderValue<Int_t> ohdu(readEvents, "ohdu");
	TTreeReaderValue<Float_t> xBary(readEvents, "xBary");
	TTreeReaderValue<Float_t> yBary(readEvents, "yBary");
	
	
	int nEvent;
	double xBaryEvent;
	double yBaryEvent;
	double xVarEvent;
	double yVarEvent;
	vector<Int_t>    xEvent;
	vector<Int_t>    yEvent;
	vector<Int_t>    maskEvent;
	vector<Float_t>  adcEvent;
	TBranch *nEventBranch = hitSumm->Branch("nEvent",&nEvent,"nEvent/I");
	TBranch *xVarEventBranch = hitSumm->Branch("xVarEvent",&xVarEvent,"xVarEvent/D");
	TBranch *yVarEventBranch = hitSumm->Branch("yVarEvent",&yVarEvent,"yVarEvent/D");
	TBranch *xBaryEventBranch = hitSumm->Branch("xBaryEvent",&xBaryEvent,"xBaryEvent/D");
	TBranch *yBaryEventBranch = hitSumm->Branch("yBaryEvent",&yBaryEvent,"yBaryEvent/D");
	TBranch *xEventBranch = hitSumm->Branch("xEvent",&xEvent[0],"xEvent[nEvent]/I");
	TBranch *yEventBranch = hitSumm->Branch("yEvent",&yEvent[0],"xEvent[nEvent]/I");
	TBranch *adcEventBranch = hitSumm->Branch("adcEvent",&adcEvent[0],"adcEvent[nEvent]/F");
	TBranch *maskEventBranch = hitSumm->Branch("maskEvent",&maskEvent[0],"maskEvent[nEvent]/I");

	int nx1e;
	while (readEvents.Next()) {
		// if (*ohdu!=hdu) {j1++; continue;}
		xEvent.clear();
		yEvent.clear();
		adcEvent.clear();
		nx1e=x1e.size();
		nEvent=0;
		
		for (int i = nx1e; i --> 0; ){
			if (*ohdu!=ohdu1e[i]) { continue;}
			for (int k = 0, nk =  xPix.GetSize(); k < nk; ++k) {
				if (pow(pow(x1e[i]-xPix[k],2)+pow(y1e[i]-yPix[k],2),0.5)<halo){
					nEvent++;
					xEvent.push_back(x1e[i]);
					yEvent.push_back(y1e[i]);
					adcEvent.push_back(ePix1e[i]);
					maskEvent.push_back(mask1e[i]);
					if (!dense){ //if dense is true, use all events, even when they repeat themselves.
						x1e.erase (x1e.begin()+ i);
						y1e.erase (y1e.begin()+ i);
						ePix1e.erase (ePix1e.begin()+ i);
					}
					break;
				}
			}
		}



		//let's compute variance of 1e- events around a track.
		hitSumm->GetEntry(j1);
		// baricenter of events around a track.
		xBaryEvent=accumulate( xEvent.begin(), xEvent.end(), 0.0) / xEvent.size();
		yBaryEvent=accumulate( yEvent.begin(), yEvent.end(), 0.0) / yEvent.size();
		
		double wSum=0;
		double xSum=0;
		double ySum=0;
		for(unsigned int i=0;i<xEvent.size();++i){
					double dx = (xEvent[i] - *xBary);
					double dy = (yEvent[i] - *yBary);
					xSum +=dx*dx*adcEvent[i];
					ySum +=dy*dy*adcEvent[i];
					wSum+=adcEvent[i];
		}

		xVarEvent = xSum/wSum;
		yVarEvent = ySum/wSum;
		hitSumm->SetBranchAddress("xEvent",&xEvent[0]);
		hitSumm->SetBranchAddress("yEvent",&yEvent[0]);
		hitSumm->SetBranchAddress("maskEvent",&maskEvent[0]);
		hitSumm->SetBranchAddress("adcEvent",&adcEvent[0]);
		xVarEventBranch->Fill();
		yVarEventBranch->Fill();
		xBaryEventBranch->Fill();
		yBaryEventBranch->Fill();
		nEventBranch->Fill();
		xEventBranch->Fill();
		maskEventBranch->Fill();
		yEventBranch->Fill();
		adcEventBranch->Fill();	
		// if (hdu==2 && nEvent>0){break;}
		j1++;
		
	}
	
	return 0;
}






int analyseHalo(const string &infile, const int &halo, const bool &dense){

	int status = 0;

	TFile *file = TFile::Open(infile.c_str());
	TTree* data = nullptr;
	file->GetObject("hitSumm",data);
	
	//set file name
	std::vector<std::string> words;split(infile,words,'/');
	std::string outfilename="output_"+words.back();
	TFile *outfile = new TFile(outfilename.c_str(),"RECREATE");


	cout << "\nProcessing file " << words.back() << " with a halo of " << halo << " pixels. Dense configuration = " << dense << ".\n"<<endl;




	//getting NCOL and NROW
	TTree* config =nullptr; file->GetObject("headerTree_0",config); config->GetEvent(0);
	
	char* NROWchar;NROWchar=((TLeafC *) config->GetBranch("NROW")->GetLeaf("string"))->GetValueString();istringstream buffer(NROWchar); buffer >> NROW; 
	char* NCOLchar;NCOLchar=((TLeafC *) config->GetBranch("NCOL")->GetLeaf("string"))->GetValueString();istringstream buffer2(NCOLchar); buffer2 >> NCOL; 
	char* CCDNROWchar;CCDNROWchar=((TLeafC *) config->GetBranch("CCDNROW")->GetLeaf("string"))->GetValueString();istringstream buffer3(CCDNROWchar); buffer3 >> CCDNROW; 
	char* CCDNCOLchar;CCDNCOLchar=((TLeafC *) config->GetBranch("CCDNCOL")->GetLeaf("string"))->GetValueString();istringstream buffer4(CCDNCOLchar); buffer4 >> CCDNCOL; 
	char* CCDNPRESchar;CCDNPRESchar=((TLeafC *) config->GetBranch("CCDNPRES")->GetLeaf("string"))->GetValueString();istringstream buffer5(CCDNPRESchar); buffer5 >> CCDNPRES; 

	//These lines will crash old files as they do not have these branches
	NBINROW=1;
	NBINCOL=1;
	if (!old){
	char* NBINROWchar;NBINROWchar=((TLeafC *) config->GetBranch("NBINROW")->GetLeaf("string"))->GetValueString();istringstream buffer6(NBINROWchar); buffer6 >> NBINROW; 
	char* NBINCOLchar;NBINCOLchar=((TLeafC *) config->GetBranch("NBINCOL")->GetLeaf("string"))->GetValueString();istringstream buffer7(NBINCOLchar); buffer7 >> NBINCOL; 
	}



	// int CCD[NCOL*NROW] = {0};
	// std::fill_n(CCD, NCOL* NROW, 0);

	//create matrix of mask
	vector<Int_t> masks(NROW*NCOL*2, 0); //2 is # of hdu's
	// if (NROW>3072){NROW=3072;}
	


	// create cuts
	std::string cutsEvents ="e>=100 && xVar>0 && yVar>0 && xMin>"+std::to_string(CCDNPRES)+" && xMax<="+std::to_string((CCDNCOL/2+CCDNPRES)/NBINCOL)+" && yMin>0 && yMax<"+std::to_string(min(NROW,CCDNROW/(2*NBINROW)))+" && ohdu<3";
	std::string cutsPix="ePix>0.63 && ePix<1.5 && !(mask & (1+4+16+128+512+1024)) && x>"+std::to_string(CCDNPRES)+" && x<="+std::to_string((CCDNCOL/2+CCDNPRES)/NBINCOL)+" && y>0 && y<"+std::to_string(min(NROW,CCDNROW/(2*NBINROW)))+" && ohdu<3";


	// Clone Tree
	TTree* hitSumm = nullptr;
	hitSumm = data->CopyTree(cutsEvents.c_str());
	hitSumm->SetName("hitSumm");
	hitSumm->SetBranchStatus("*",1);
	// create TTreeReader for that Tree
	TTreeReader readEvents("hitSumm", outfile);



	

	

	// Select 1e- events
	vector<Int_t>    x1e;
	vector<Int_t>    y1e;
	vector<Double_t>    ePix1e;
	vector<Int_t>    ohdu1e;
	vector<Int_t>    mask1e;
	// cutsPix.pop_back();
	// char charhdu = hdu+48; //cast int to char
	// cutsPix.push_back(charhdu);
	cout << "Cut for 1e events = " << cutsPix <<"\n"<<endl;
	Compute1eEvents(file, outfile, x1e, y1e, ePix1e, ohdu1e,mask1e, cutsPix,masks);

	// Select events (tracks)
	cout << "Cut for tracks = " << cutsEvents << "\n"<<endl;
	cutEvents(file, outfile, readEvents, hitSumm, masks);


	// writes tree with new branches
	WriteTree(readEvents, hitSumm, x1e, y1e, ePix1e,ohdu1e,mask1e);
	x1e.clear();
	y1e.clear();
	ePix1e.clear();
	
	
	

	//Clone config tree (skViewer needs it)
	TTree* imgParTreeConfig = nullptr;
	file->GetObject("imgParTree",imgParTreeConfig);
	TTree* imgParTree = nullptr;
	imgParTree = imgParTreeConfig->CloneTree();
	imgParTree->SetName("imgParTree");
	imgParTree->SetBranchStatus("*",1);
	delete imgParTreeConfig;

	
	//add halo size to imgParTree
	int haloVal;
	TBranch *haloBranch = imgParTree->Branch("halo",&haloVal,"halo/I");
	haloVal=halo;
	haloBranch->Fill();
	
	outfile->Write();
	outfile->Close();
	
	return status;
}




/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[]){


time_t start,end;
double dif;
time (&start);

std::string infile(argv[1]);
dense=false; //global variable
old=false; //global variable
halo=30; //global variable

int returnCode = processCommandLineArgs( argc, argv, dense, halo); //processing arguments

int status = analyseHalo(infile, halo, dense);


time (&end);
dif = difftime (end,start);

cout << "\n\n---------------------------------" <<endl;
cout << "\n\nTook me " << dif << " seconds.\n\n" << endl;
cout << "---------------------------------\n\n" <<endl;
return 0;

}    // end
