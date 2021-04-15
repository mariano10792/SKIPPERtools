/* Version 28/02/2021
Autor: Dario Rodrigues
mail: dariorodriguesfm@gmail.com

Código escrito para simular mediciones de la letra LAMBDA por Microscopia Quantica con Skipper-CCD

- Coordenada (xx,yy) donde pueden llegar fotones
- Eventos de uno, dos o mas electrones generados por ruido de lectura.
- Distribucion de los pixeles que conforman cada cluster

para correrlo escribir por ejemplo ./CCD_simulation.exe 1296 0.5 0.05 1
donde el primer numero representa el numero de fotones a simular
el segundo es el ruido en electrones por pixel que hubo durante la medicion
el tercero es la fraccion de Stray Light sobre el numero total (Stray mas simulados)
el ultimo es simplemente un numero para etiquetar el RUN y poder luego diferenciar los outputs files.

Importante!!! si no se cambia el numero de RUN las nuevas salidas pisan las anteriores
*/

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

// Main control variables
int nsim=1; //number of simulations
int factor=1; //multiplicative factor
bool simulation=true; //To set calPixTree1e as a simulation, too.
// 1e cuts. Very low contamination.
double lowercut=.85; //1.5 for 2e-
double uppercut=1.5; //2.5 for 2e-
// set cuts and edgedist alias
string cuts="!(mask & (1+16+128+256+512+1024+4096+16384)) && min(bleedX,bleedY)>=100 && min(distance,edgedist)>="; //"!(mask & (1+4+16+128+256+512+1024+4096)) && distance1e>=20 && min(distance,edgedist)>=" for 2e-


// Masfile configuration variables // 
int NROW;int NCOL;int NBINROW;int NBINCOL;int CCDNCOL;int CCDNROW;int CCDNPRES;int RUNID;int LTANAME;


//////////////////Read maskfile configuration ////////////////////////////

void readConfig(TFile *file)
{
	TTree* config =nullptr; file->GetObject("headerTree_0",config); config->GetEvent(0);
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
	char* LTANAMEchar;LTANAMEchar=((TLeafC *) config->GetBranch("LTANAME")->GetLeaf("string"))->GetValueString();istringstream buffer9(LTANAMEchar); buffer9 >> LTANAME; 
}

////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[]){

	TStopwatch t; // time counter
	t.Start();

	std::vector<string> inFileList;
	for(int i=optind; i<argc-1; ++i){ // -1 to delete option -h or -p. Not final solution.
		inFileList.push_back(argv[i]);
	}

	//Initialize simulation
	int nsimextra=0; //to write a simulation into calPixTree1e
	if (simulation==true)
	{
		nsimextra++;
	}
	vector<TTree> trees(nsim+nsimextra); //output simulated trees
	TList * list = new TList; //ouptut real trees

	//Set output file name
	//set file name
	std::vector<std::string> words;
	split(inFileList[0],words,'/');
	std::string outfilename="sim_halo";
	if (simulation==true)
	{
		outfilename="sim_vs_sim_halo";
	}
	outfilename+=argv[argc-1];
	outfilename+="_";
	outfilename+=words.back();
	TFile *outfile = new TFile(("files/"+outfilename).c_str(),"RECREATE"); //output file
	
	
	
	
	for (size_t iFile = 0; iFile < inFileList.size(); iFile++)
	{
		
		// extract TTree from TFiles and read configuration
		TFile *file = TFile::Open(inFileList[iFile].c_str());
		cout << "Maskfile: " << inFileList[iFile] << endl;
		TTree *calPixTree = (TTree*)file->Get("calPixTree");
		readConfig(file);
		cuts+=argv[argc-1]; //last arg is halo
		calPixTree->SetAlias("edgedist",("min(min(x-"+std::to_string(CCDNPRES+1)+", "+std::to_string((CCDNCOL/2+CCDNPRES)/NBINCOL)+"-x), min(y-1, "+std::to_string(min(NROW,CCDNROW/(2*NBINROW)))+"-y))").c_str());


		

		// extract real 0 and 1e- events into a TTreeReader
		outfile->cd();
		TTree *calPixTreeAux = calPixTree->CopyTree((cuts+"&& ePix<"+std::to_string(uppercut)).c_str());
		calPixTreeAux->SetName("calPixTreeAux");
		calPixTreeAux->SetBranchStatus("*",1);
		TTreeReader readValues("calPixTreeAux", outfile);
		TTreeReaderValue<Int_t> ohduReader(readValues, "ohdu");
		TTreeReaderValue<Int_t> xReader(readValues, "x");
		TTreeReaderValue<Int_t> yReader(readValues, "y");
		TTreeReaderValue<Double_t> ePixReader(readValues, "ePix");




		// create xy vector with all data
		vector<vector<Int_t>> xy(5); // ohdu, x position, y position, x+y*NCOL position and boolean (1 if occupied 0 if not)
		vector <Int_t> events(4,0); // # of events per ohdu
		// extract events into xy
		int prevohdu=1;
		int size=0;
		while (readValues.Next()) {
			xy[0].push_back(*ohduReader);
			size++;
			if (prevohdu!=*ohduReader){events[prevohdu-1]+=size; prevohdu=*ohduReader; size=0;}
			xy[1].push_back(*xReader);
			xy[2].push_back(*yReader);
			if (*ePixReader>lowercut && *ePixReader<uppercut)
			{
				xy[3].push_back(1);
				xy[4].push_back(*xReader+*yReader*NCOL);
			}else{
				xy[3].push_back(0);
				xy[4].push_back(*xReader+*yReader*NCOL);
			}
		}	
		events[prevohdu-1]=size;


		//extract real x events into a TList
		if (simulation==false)
		{
			list->Add(calPixTree->CopyTree((cuts+" && ePix>"+std::to_string(lowercut)+" && ePix<"+std::to_string(uppercut)).c_str()));
		}
		
		int initcounter=0;

		// vector of TF1 we will use to fit each kolmogorov's cdf from each quad
		vector <TF1> f1(4);
		TCanvas * canvas   = new TCanvas("Fit","Fit",0,0,1200,800); //coord top lef, size x, sizey
		// canvas->cd();
		canvas->Print(("files/LadderFit_"+outfilename+".pdf[").c_str());
		
		for (size_t hdu = 1; hdu < 5; hdu++)
		{
			//create kolmogorov cdf
			Double_t *ladder = new Double_t[events[hdu-1]];
			Double_t *ladderX = new Double_t[events[hdu-1]];
			int counter=0;
			for (size_t i = 0; i < events[hdu-1]; i++)
			{
				if (xy[3][i+initcounter]==1){counter++;}
				ladder[i]=float(counter);
				ladderX[i]=i;
			}
			initcounter+=events[hdu-1];
			cout << "I counted " << counter << " events in ohdu " << hdu << "." <<  endl;

			//fit kolmogorov cdf
						
			TGraphErrors *graph = new TGraphErrors(events[hdu-1], ladderX, ladder);
			TF1 *a = new TF1(("f"+std::to_string(hdu)).c_str(), "[0]*x*x+[1]*x", 0, events[hdu-1]);
			graph->Draw();
			graph->Fit(("f"+std::to_string(hdu)).c_str(),"Q");
			f1[hdu-1] = (*a);
			// canvas->SaveAs(("LadderFit"+std::to_string(hdu)+"_"+outfilename+".pdf").c_str());
			canvas->cd();
			canvas->Print(("files/LadderFit_"+outfilename+".pdf").c_str());
			canvas->Clear();
			delete[] ladder; delete[] ladderX;
		}
		canvas->cd();
		canvas->Print(("files/LadderFit_"+outfilename+".pdf]").c_str());
		
		//simulate readout DC and save to trees
		Int_t N1; Int_t x; Int_t y; Int_t ohdu;
		TRandom3 DC(0);
		
		for (size_t t = 0; t < nsim; t++)
		{
			trees[t].SetName(("simPixTree"+std::to_string(t)).c_str());
			trees[t].Branch("RUNID",&t,"RUNID/I");
			trees[t].Branch("ohdu",&ohdu,"ohdu/I");
			trees[t].Branch("x",&x,"x/I");
			trees[t].Branch("y",&y,"y/I");

			for (size_t hdu = 1; hdu < 5; hdu++)
			{
				ohdu=hdu;
				for (size_t i = 0; i < events[hdu-1]; i++)
				{
					if (DC.Poisson((f1[hdu-1].Eval(i)-f1[hdu-1].Eval(i-1))*factor)>0)
					{
						x=xy[1][i];
						y=xy[2][i];
						trees[t].Fill();
					}
				}
			}
			
		}
		if (simulation==true)
		{
			int x; int y; int ohdu;
			trees[nsim].SetName("calPixTree1e");
			trees[nsim].Branch("RUNID",&nsim,"RUNID/I");
			trees[nsim].Branch("ohdu",&ohdu,"ohdu/I");
			trees[nsim].Branch("x",&x,"x/I");
			trees[nsim].Branch("y",&y,"y/I");

			for (size_t hdu = 1; hdu < 5; hdu++)
			{
				ohdu=hdu;
				for (size_t i = 0; i < events[hdu-1]; i++)
				{
					if (DC.Poisson((f1[hdu-1].Eval(i)-f1[hdu-1].Eval(i-1))*factor)>0)
					{
						x=xy[1][i];
						y=xy[2][i];
						trees[nsim].Fill();
					}
				}
			}

		}

		delete calPixTreeAux; // delete so future calPixTreeAux is re-used
	}

	//simulation header
	TTree * simHeader = new TTree("simHeader","simHeader");
	simHeader->Branch("nsim",&nsim,"nsim/I");
	simHeader->Branch("LTANAME",&LTANAME,"LTANAME/I");
	simHeader->Branch("runID",&RUNID,"runID/I");
	simHeader->Branch("factor",&factor,"factor/I");
	simHeader->Fill();


	//after filling write each tree
	for (size_t t = 0; t < nsim+nsimextra; t++)
	{trees[t].Write();}


	// merge TList with real x events or write simuated tree into calPixTree1e
	if (simulation==false)
	{
		TTree *calPixTree1e = TTree::MergeTrees(list);
		calPixTree1e->SetName("calPixTree1e");
		calPixTree1e->Write();
	}
	
	
	// Write output file
	outfile->Write();
	outfile->Close();

t.Stop();
t.Print();

return 0;
}    // end program