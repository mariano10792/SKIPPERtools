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
#include <fstream>
#include <chrono>
#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TStopwatch.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TLeafC.h"
#include "TStyle.h"
#include "TRandom3.h"
#include "TChain.h"
#include <string>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <numeric>

#include "gConfig.h"
#include "gCal.h"
#include "tinyxml2.h"
#include <unistd.h>

using namespace std;
long totpix;
int run = 0;
Int_t xSize     = -1;
Int_t ySize     = -1;
gConfig &gc = gConfig::getInstance();
gCal &gcal = gCal::getInstance();


////////////////////////////////////////////////////////
// Function used to read fits files

std::string trim(const std::string& str, const std::string& whitespace = " \t\'"){ // removes leading and trailing spaces
    const auto strBegin = str.find_first_not_of(whitespace);
    if (strBegin == std::string::npos) return ""; // no content
    const auto strEnd = str.find_last_not_of(whitespace);
    const auto strRange = strEnd - strBegin + 1;
    return str.substr(strBegin, strRange);
}

int getExpoInfoFromHeader(fitsfile *fptr, Int_t &tShut, Int_t &tExpo){
    int status = 0;   /*  CFITSIO status value MUST be initialized to zero!  */

    int hdutype;
    char keyValue[1024] = "";
    char comment[1024]  = "";

    Int_t tDate = -1;
    tShut = -1;
    tExpo = -1;
    // Get general data from ext=0 (hdu=1) header
    fits_movabs_hdu(fptr, 1, &hdutype, &status);
    if(status!=0) return -2;

    fits_read_keyword(fptr, "DATE", keyValue, comment, &status);
    if(status==0){ // key exist
        std::tm tm = {};
        strptime(trim(keyValue).c_str(), "%Y-%m-%dT%H:%M:%S", &tm);
        time_t t = std::mktime(&tm);
        tDate = t;	    
    }

    fits_read_keyword(fptr, "UTSHUT", keyValue, comment, &status);
    if(status==0){ // key exist
        std::tm tm = {};
        strptime(trim(keyValue).c_str(), "%Y-%m-%dT%H:%M:%S", &tm);
        time_t t = std::mktime(&tm);
        tShut = t;	    
    }
    tExpo = tDate - tShut;

    return 0;
}

bool readCardValue(fitsfile  *infptr, const char *keyName, double &value){

    int status = 0;
    char record[1024] = "";
    fits_read_card(infptr, keyName, record, &status);
    if(status==KEY_NO_EXIST){
        status=0;
        return false;
    }
    else{
        string sRec(record);
        size_t tPosEq = sRec.find("=");
        size_t tPosSl = sRec.find("/");
        string sVal(sRec.substr(tPosEq+1, tPosSl-tPosEq-1));
        std::replace( sVal.begin(), sVal.end(), '\'', ' ');
        istringstream recISS( sVal );
        recISS >> value;
        return true;
    }

}

bool readIntCard(fitsfile  *infptr, const char *keyName, int &value){
    int status = 0;
    int oldValue = value;
    fits_read_key(infptr, TINT, keyName, &value, NULL, &status);
    if (status!=0) value = oldValue; //if any error (no key, or value empty/invalid), restore the original value
    return (status!=KEY_NO_EXIST);
}
////////////////////////////////////////////////////////

bool alone(int x0, double* &dummyimageArray)
{
    bool status=0;
    if (dummyimageArray[x0]==0 && dummyimageArray[x0+1]==0 && dummyimageArray[x0-1]==0 && dummyimageArray[x0+xSize]==0 && dummyimageArray[x0+xSize+1]==0 && dummyimageArray[x0+xSize-1]==0 && dummyimageArray[x0-xSize]==0 && dummyimageArray[x0-xSize-1]==0 && dummyimageArray[x0-xSize+1]==0 && dummyimageArray[x0+2]==0 && dummyimageArray[x0-2]==0 && dummyimageArray[x0-2*xSize]==0 && dummyimageArray[x0+2*xSize]==0 && dummyimageArray[x0+2*xSize+1]==0 && dummyimageArray[x0+2*xSize+2]==0 && dummyimageArray[x0+2*xSize-2]==0 && dummyimageArray[x0+2*xSize-1]==0 && dummyimageArray[x0-2*xSize+1]==0 && dummyimageArray[x0-2*xSize+2]==0 && dummyimageArray[x0-2*xSize-2]==0 && dummyimageArray[x0-2*xSize-1]==0 && dummyimageArray[x0+xSize+2]==0 && dummyimageArray[x0+xSize-2]==0 && dummyimageArray[x0-xSize+2]==0 && dummyimageArray[x0-xSize-2]==0)
    {}
    else
    {
        status=1;
    }
    
    return status;
}

int around(int x0, double* &dummyimageArray)
{
    if (dummyimageArray[x0]==0){dummyimageArray[x0]+=2;}
    if (dummyimageArray[x0+1]==0){dummyimageArray[x0+1]+=2;}
    if (dummyimageArray[x0-1]==0){dummyimageArray[x0-1]+=2;}
    if (dummyimageArray[x0+1+xSize]==0){dummyimageArray[x0+1+xSize]+=2;}
    if (dummyimageArray[x0-1+xSize]==0){dummyimageArray[x0-1+xSize]+=2;}
    if (dummyimageArray[x0+xSize]==0){dummyimageArray[x0+xSize]+=2;}
    if (dummyimageArray[x0-1-xSize]==0){dummyimageArray[x0-1-xSize]+=2;}
    if (dummyimageArray[x0+1-xSize]==0){dummyimageArray[x0+1-xSize]+=2;}
    if (dummyimageArray[x0-xSize]==0){dummyimageArray[x0-xSize]+=2;}
    return 0;    
}

int copyStructure(const char* inF, const char *outF){

    fitsfile  *outfptr; /* FITS file pointers defined in fitsio.h */
    fitsfile *infptr;   /* FITS file pointers defined in fitsio.h */

    int status = 0;
    // int single = 0;

    int hdutype, bitpix, naxis = 0, nkeys;
    int nhdu = 0;
    long naxes[9] = {1, 1, 1, 1, 1, 1, 1, 1, 1};
    long totpix = 0;

    fits_open_file(&infptr, inF, READONLY, &status); /* Open the input file */
    if (status != 0) return(status);

    fits_get_num_hdus(infptr, &nhdu, &status);  

    fits_create_file(&outfptr, outF, &status);/* Create the output file */
    if (status != 0) return(status);


    for (int n=1; n<=nhdu; ++n)  /* Main loop through each extension */
    { 
        /* get image dimensions and total number of pixels in image */
        fits_movabs_hdu(infptr, n, &hdutype, &status);
        for (int i = 0; i < 9; ++i) naxes[i] = 1;
        fits_get_img_param(infptr, 9, &bitpix, &naxis, naxes, &status);
        totpix = naxes[0] * naxes[1] * naxes[2] * naxes[3] * naxes[4] * naxes[5] * naxes[6] * naxes[7] * naxes[8];

        if (hdutype != IMAGE_HDU || naxis == 0 || totpix == 0){
            /* just copy tables and null images */
            fits_copy_hdu(infptr, outfptr, 0, &status);
            if (status != 0) return(status);
        }
        else{
            fits_create_img(outfptr, -32, naxis, naxes, &status);//image of signed int values
            if (status != 0) return(status);

            /* copy the header keywords */
            fits_get_hdrspace(infptr, &nkeys, NULL, &status); 
            for (int i = 1; i <= nkeys; ++i){
                char card[FLEN_CARD];
                fits_read_record(infptr, i, card, &status);
                if (fits_get_keyclass(card) > TYP_CMPRS_KEY) fits_write_record(outfptr, card, &status);
            }

        }
    }

    fits_close_file(infptr, &status);
    fits_close_file(outfptr,  &status);

    return status;
}

bool fileExist(const char *fileName){
    ifstream in(fileName,ios::in);

    if(in.fail()){
        //cout <<"\nError reading file: " << fileName <<"\nThe file doesn't exist!\n\n";
        in.close();
        return false;
    }

    in.close();
    return true;
}

int deleteFile(const char *fileName){
    cout << "Will overwrite: " << fileName << endl << endl;
    return unlink(fileName);
}

int processCommandLineArgs(const int argc, char *argv[], string &maskFile, vector<string> &inFileList){

    if(argc == 1) return 1;
    int opt=0;
    while ( (opt = getopt(argc, argv, "c:C:m:n:q?")) != -1) {
        switch (opt) {
            case 'm':
                maskFile = optarg;
                break;
            case 'n':
                run = atoi(optarg);
                break;
            case 'C':
                if(gcal.readConfFile(optarg) == false){
                    return 1;
                }
                break;
            case 'c':
                if(gc.readConfFile(optarg) == false){
                    return 1;
                }
                break;
            case 'q':
                break;
            default: /* '?' */
                return 1;
        }
    }

    for(int i=optind; i<argc; ++i){
        inFileList.push_back(argv[i]);
        if(!fileExist(argv[i])){
            cout << "\nError reading input file: " << argv[i] <<"\nThe file doesn't exist!\n\n";
            return 1;
        }
    }

    if(inFileList.size()==0){
        cerr << "Error: no input file(s) provided!\n\n";
        return 1;
    }

    return 0;
}

int readMask(string maskFile, vector <int*> &masks){
    int status = 0;
    int nhdu = 0;
    double nulval = 0.;
    int anynul = 0;
    
    const char* maskName = maskFile.c_str();

    fitsfile  *infptr; /* FITS file pointers defined in fitsio.h */
    fits_open_file(&infptr, maskName, READONLY, &status); /* Open the input file */
    if (status != 0) return(status);
    fits_get_num_hdus(infptr, &nhdu, &status);
    if (status != 0) return(status);


    for (unsigned int eN=1; eN<=nhdu; ++eN)  /* Main loop through each extension */
    {
        const int n = eN;

        /* get input image dimensions and total number of pixels in image */
        int hdutype, bitpix, naxis = 0;
        long naxes[9] = {1, 1, 1, 1, 1, 1, 1, 1, 1};
        fits_movabs_hdu(infptr, n, &hdutype, &status);
        for (int i = 0; i < 9; ++i) naxes[i] = 1;
        fits_get_img_param(infptr, 9, &bitpix, &naxis, naxes, &status);
        totpix = naxes[0] * naxes[1];

        /* Don't try to process data if the hdu is empty */    
        if (hdutype != IMAGE_HDU || naxis == 0 || totpix == 0){
            masks.push_back(0);
            continue;
        }

        int* maskArray = new int[totpix];

        /* Open the input file */
        fits_movabs_hdu(infptr, n, &hdutype, &status);
        if (status != 0) return(status);

        /* Read the images as doubles, regardless of actual datatype. */
        long fpixel[2]={1,1};
        long lpixel[2]={naxes[0],naxes[1]};
        long inc[2]={1,1};
        fits_read_subset(infptr, TINT, fpixel, lpixel, inc, &nulval, maskArray, &anynul, &status);
        if (status != 0){
            fits_report_error(stderr, status);
            return(status);
        }
        masks.push_back(maskArray);
    }
    
    fits_close_file(infptr, &status);

    return status;
}


int writeImage(const char *maskfName, const int ext, const long totpix, double* imageArray,int counter) {
    int status = 0;
    /* Open the output file */
    fitsfile  *maskfptr; /* FITS file pointers defined in fitsio.h */
    fits_open_file(&maskfptr, maskfName, READWRITE, &status);
    if (status != 0) return(status);

    // fits_update_key(maskfptr, TINT, "EVENTS", &counter,"EVENTS", &status);
    fits_update_key(maskfptr, TINT, "CHID", &counter,"CHID", &status);

    int hdutype;
    fits_movabs_hdu(maskfptr, ext, &hdutype, &status);
    long fpixel[2]={1,1};
    fits_write_pix(maskfptr, TDOUBLE, fpixel, totpix, imageArray, &status);
    fits_close_file(maskfptr,  &status);

    return status;
}


int extractImages(const vector<string> &inFileList, string maskFile, vector <int*> &masks, vector <double*> &inputimages){
    int status = 0;
    double nulval = 0.;
    int anynul = 0;

    // read the provided mask file, if there is one
    readMask(maskFile, masks);
    // End of mask handling

    int nhdu = 0;
    const unsigned int nFiles  = inFileList.size();

    // Card values to be read from fits file
    Int_t tShut     = -1;
    Int_t tExpo     = -1;
    Int_t tRunID    = -1;
    Int_t tLTANAME    = -1;
    Double_t tTemp  = -1;


    // Input image handling
    
    for(unsigned int fn=0; fn < nFiles; ++fn){

        fitsfile  *infptr; /* FITS file pointers defined in fitsio.h */
        fits_open_file(&infptr, inFileList[fn].c_str(), READONLY, &status); /* Open the input file */
        if (status != 0) return(status);
        fits_get_num_hdus(infptr, &nhdu, &status);
        if (status != 0) return(status);

        getExpoInfoFromHeader(infptr, tShut, tExpo);
        readIntCard(infptr, "RUNID", tRunID);
        readIntCard(infptr, "LTANAME", tLTANAME);
        readCardValue(infptr, "TEMPER", tTemp);

        for (unsigned int eN=1; eN<=nhdu; ++eN)  /* Main loop through each extension */
        {

            const int n = eN;

            /* get input image dimensions and total number of pixels in image */
            int hdutype, bitpix, naxis = 0;
            long naxes[9] = {1, 1, 1, 1, 1, 1, 1, 1, 1};
            fits_movabs_hdu(infptr, n, &hdutype, &status);
            for (int i = 0; i < 9; ++i) naxes[i] = 1;
            fits_get_img_param(infptr, 9, &bitpix, &naxis, naxes, &status);
            long totpix = naxes[0] * naxes[1];

            /* Don't try to process data if the hdu is empty */    
            if (hdutype != IMAGE_HDU || naxis == 0 || totpix == 0){
                continue;
            }

            double* ePixArray = new double[totpix]; //electrons per pixel, not rounded

            /* Open the input file */
            fits_movabs_hdu(infptr, n, &hdutype, &status);
            if (status != 0) return(status);
            if(xSize<naxes[0]) xSize = naxes[0];
            if(ySize<naxes[1]) ySize = naxes[1];

            /* Read the images as doubles, regardless of actual datatype. */
            long fpixel[2]={1,1};
            long lpixel[2]={naxes[0],naxes[1]};
            long inc[2]={1,1};

            // read the raw ADUs into ePixArray
            fits_read_subset(infptr, TDOUBLE, fpixel, lpixel, inc, &nulval, ePixArray, &anynul, &status);
            if (status != 0) return(status);

            // const double kCal = gcal.isValid()?gcal.getExtCal(eN):gc.getExtCal(eN);
            // for (int i=0; i<totpix; ++i) { //scale ePixArray to get electrons
            //     ePixArray[i] /= kCal;
            // }
            // cout << kCal << endl;

            cout << "\nProcessing runID " << tRunID << " ohdu " << eN << ":\n";

            inputimages.push_back(ePixArray);
        }

        /* Close the input file */
        fits_close_file(infptr,  &status);   

    }

    return status;
}


int main(int argc, char *argv[])
{
    time_t start,end;
    double dif;
    time (&start);

    vector<string> inFileList;
    string maskFile;
    int returnCode = processCommandLineArgs(argc,argv,maskFile,inFileList);
    if(returnCode!=0){
        return returnCode;
    } 

    vector <int*> masks;
    vector <double*> inputimages;

    string simfName="";
    const int nameStart     = std::max((int)(inFileList[0].rfind("/")+1), 0);
    const int nameEnd       = std::max((int)(inFileList[0].rfind(".fits")+1), 0);
    const int extraProcName = (inFileList[0].substr(nameStart,5)=="proc_") ? 5 : 0;
    string outFileBaseName  = inFileList[0].substr(nameStart+extraProcName, nameEnd-nameStart-extraProcName-1);
    simfName = "sim_"+outFileBaseName+"_"+std::to_string(run)+".fits";
    if(fileExist(simfName.c_str())) deleteFile(simfName.c_str());
    copyStructure(inFileList[0].c_str(), simfName.c_str());

    extractImages(inFileList, maskFile, masks, inputimages);

    vector <double*> outputimages;
    vector <int*> dummyimages;

    TRandom3 rndm(0);
    int x0;
    for (Int_t i = 0; i < masks.size(); i++)
    {
        const double kCal = gcal.isValid()?gcal.getExtCal(i+1):gc.getExtCal(i+1);
        double* outputimageArray = new double[totpix];
        std::fill_n(outputimageArray, totpix, 0.0);
        double* dummyimageArray = new double[totpix];
        std::fill_n(dummyimageArray, totpix, 1);
        int* maskArray = masks[i];
        double* inputimageArray = inputimages[i];

        // write image
        for (size_t j = 0; j< totpix; j++)
        {
            if (maskArray[j]==0)
            {
                outputimageArray[j]=inputimageArray[j];
                dummyimageArray[j]=0;
            }
        }
      
        int counter=0;
        int notcounter=0;
        while (true)
        {
            x0 = round(rndm.Uniform(xSize,totpix-xSize));
            if (alone(x0,dummyimageArray) && alone(x0-xSize+1,dummyimageArray)==0)
            {
                outputimageArray[x0]+=kCal;
                outputimageArray[x0-xSize+1]+=kCal;
                around(x0,dummyimageArray);
                around(x0-xSize+1,dummyimageArray);
                counter++;
                notcounter=0;
            }
            else if (dummyimageArray[x0]==2 || dummyimageArray[x0-xSize+1]==2) notcounter++;
            if (notcounter>10000) break;//si no pudiste posicionar un evento 10000 veces seguidas, desistí.
        }

        writeImage(simfName.c_str(),i+1,totpix,outputimageArray,counter);
        delete[] outputimageArray;
        delete[] inputimageArray;
        delete[] dummyimageArray;
        delete[] maskArray;
    }
    

    time (&end);
    dif = difftime (end,start);
    cout << "\nAll done!\n" << "-> It took me " << dif << " seconds to do it!\n\n";


    return 0;
}

//make && ./DM.exe /mnt/mariano/MINOSnewfirmware/commissioning_data/minos3/proc_corr_proc_skp_moduleC40_41-ssc16_17-lta20_60_TEMP135K-run1_NROW520_NBINROW1_NCOL3200_NBINCOL1_EXPOSURE72000_CLEAR1800_1_250.fits -m /mnt/mariano/MINOSnewfirmware/commissioning_data/minos3/mask_corr_proc_skp_moduleC40_41-ssc16_17-lta20_60_TEMP135K-run1_NROW520_NBINROW1_NCOL3200_NBINCOL1_EXPOSURE72000_CLEAR1800_1_250.fits -c extractConfig.xml -C /mnt/mariano/MINOSnewfirmware/commissioning_data/minos3/cal_proc_corr_proc_skp_moduleC40_41-ssc16_17-lta20_60_TEMP135K-run1_NROW520_NBINROW1_NCOL3200_NBINCOL1_EXPOSURE72000_CLEAR1800_1_250.xml -n 0 && ds9 sim_corr_proc_skp_moduleC40_41-ssc16_17-lta20_60_TEMP135K-run1_NROW520_NBINROW1_NCOL3200_NBINCOL1_EXPOSURE72000_CLEAR1800_1_250_0.fits