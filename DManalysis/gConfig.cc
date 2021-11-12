#include "tinyxml2.h"
#include "gConfig.h"
#include "globalConstants.h"
#include <iostream>
#include <sstream>


using namespace std;

gConfig::gConfig()
    : fDefaultCal(-1)
    , fEpixCut(0.5)
    , fBleedXCut(50)
    , fBleedYCut(50)
    , fHaloCut(50)
    , fEdgeCut(-1)
    , fClusterThr(5)
    , fClusterCut(4)
    , fClusterIgnore(-1)
    , fCrosstalkThr(700)
    , fHasCal(false)
    , fIsValid(false)
    , fLooseCluster(false)
    , fLooseClustersNoNeighbors(false)
    , fUseDiagonalPix(false)
    , fMaskSmall(false)
      , fSaveTracks(false)
{}

float gConfig::getExtCal(const int ext){

    if (fExt2Cal.count(ext)==0) {
        return fDefaultCal;
    }
    return fExt2Cal[ext];
}

bool gConfig::isBadPix(const int ohdu, const int x, const int y){
    if (fExt2BadPix.count(ohdu)==0) {
        return false;
    } else {
        return fExt2BadPix[ohdu].count(make_pair(x,y))!=0;
    }
}

bool gConfig::isBadCol(const int ohdu, const int x){
    if (fExt2BadCol.count(ohdu)==0) {
        return false;
    } else {
        return fExt2BadCol[ohdu].count(x)!=0;
    }
}

bool gConfig::isBleedCol(const int ohdu, const int x){
    if (fExt2BleedCol.count(ohdu)==0) {
        return false;
    } else {
        return fExt2BleedCol[ohdu].count(x)!=0;
    }
}

bool gConfig::isPastBleedXEdge(const int ohdu, const int x){
    if (fExt2BleedXEdge.count(ohdu)==0) {
        return false;
    } else {
        return x>=fExt2BleedXEdge[ohdu];
    }
}

void gConfig::printVariables(){
    cout << endl;
    cout << "===== Configuration parameters =====\n";
    cout << "badPixels:\n";
    typedef std::map<int, std::set<std::pair<int,int>>>::iterator it_badpix;
    for(it_badpix it = fExt2BadPix.begin(); it != fExt2BadPix.end(); it++) {
        typedef std::set<std::pair<int,int>>::iterator it_type;
        for(it_type it2 = it->second.begin(); it2 != it->second.end(); it2++) {
            cout << "\text" << it->first << " x: " << it2->first << " y: " << it2->second << endl;
        }
    }

    typedef std::map<int, std::set<int>>::iterator it_col;

    cout << "badCols:\n";
    for(it_col it = fExt2BadCol.begin(); it != fExt2BadCol.end(); it++) {
        typedef std::set<int>::iterator it_type;
        for(it_type it2 = it->second.begin(); it2 != it->second.end(); it2++) {
            cout << "\text" << it->first << " x: " << *it2 << endl;
        }
    }

    cout << "bleedCols:\n";
    for(it_col it = fExt2BleedCol.begin(); it != fExt2BleedCol.end(); it++) {
        typedef std::set<int>::iterator it_type;
        for(it_type it2 = it->second.begin(); it2 != it->second.end(); it2++) {
            cout << "\text" << it->first << " x: " << *it2 << endl;
        }
    }

    cout << "bleedXEdges:\n";
    for(std::map<int, int>::iterator it = fExt2BleedXEdge.begin(); it != fExt2BleedXEdge.end(); it++) {
        cout << "\text" << it->first << " x: " << it->second << endl;
    }

    cout << "cuts:\n";
    cout << "\tePix:          " << fEpixCut << endl;
    cout << "\tbleedX:        " << fBleedXCut << endl;
    cout << "\tbleedY:        " << fBleedYCut << endl;
    cout << "\thalo:          " << fHaloCut << endl;
    cout << "\tedgeCut:       " << fEdgeCut << endl;
    cout << "\tclusterThr:    " << fClusterThr << endl;
    cout << "\tclusterCut:    " << fClusterCut << endl;
    cout << "\tclusterIgnore: " << fClusterIgnore << endl;
    cout << "\tcrosstalkThr:  " << fCrosstalkThr << endl;

    cout << "\tlooseCluster:              " << fLooseCluster             << endl;
    cout << "\tlooseClustersNoNeighbors:  " << fLooseClustersNoNeighbors << endl;
    cout << "\tuseDiagonalPix:            " << fUseDiagonalPix           << endl;
    cout << "\tmaskSmall:                 " << fMaskSmall                << endl;  

    cout << "extra:\n";
    cout << "\tsaveTracks:                " << fSaveTracks               << endl;
    cout << "\tsaveTrackCuts:             " << fTracksCuts               << endl;  

    cout << "====================================\n";
}

bool gConfig::readConfFile(const char* confFileName = "extractConfig.xml"){
    fFilename = confFileName;

    tinyxml2::XMLDocument doc;
    if(doc.LoadFile( confFileName ) != 0){
        cerr << red << "\nCan't read config file! Will not continue.\n\n" << normal;
        return false;
    }
    //   doc.Print();

    tinyxml2::XMLElement* config = doc.FirstChildElement("extractConfig");
    if ( config == NULL ){
        cerr << "Missing \'extractConfig\' in config file.\n\n";
        return false;
    }

    /* Calibration */
    tinyxml2::XMLElement* calibration = config->FirstChildElement("calibration");
    if( calibration != NULL ){ //it's OK if this is missing (the calibration should be specified with -C option)
        tinyxml2::XMLElement* hdu = calibration->FirstChildElement("hdu");
        while (hdu != NULL) {
            int hdunum;
            if (hdu->QueryIntAttribute("num", &hdunum) != tinyxml2::XML_SUCCESS) {
                cerr << "Missing or invalid \'num\' attribute in config file.\n\n";
                return false;
            }
            float gain;
            if (hdu->QueryFloatAttribute("gain", &gain) != tinyxml2::XML_SUCCESS) {
                cerr << "Missing or invalid \'gain\' attribute in config file.\n\n";
                return false;
            }
            fExt2Cal[hdunum] = gain;
            hdu = hdu->NextSiblingElement("hdu");
        }

        tinyxml2::XMLElement* defaultCal = calibration->FirstChildElement("default");
        if( defaultCal == NULL ){
            cerr << "Missing \'default\' in config file.\n\n";
            return false;
        }
        if (defaultCal->QueryFloatAttribute("gain", &fDefaultCal) != tinyxml2::XML_SUCCESS) {
            cerr << "Missing or invalid \'gain\' attribute in config file.\n\n";
            return false;
        }
        fHasCal = true;
    }

    /* badPixels */
    tinyxml2::XMLElement* badPixels = config->FirstChildElement("badPixels");
    if( badPixels == NULL ){
        cerr << "Missing \'badPixels\' in config file.\n\n";
    } else {
        tinyxml2::XMLElement* pixel = badPixels->FirstChildElement("pixel");
        while (pixel != NULL) {
            int hdu,x,y;
            if (pixel->QueryIntAttribute("hdu", &hdu) != tinyxml2::XML_SUCCESS) {
                cerr << "Missing or invalid \'hdu\' attribute in config file.\n\n";
                return false;
            }
            if (pixel->QueryIntAttribute("x", &x) != tinyxml2::XML_SUCCESS) {
                cerr << "Missing or invalid \'x\' attribute in config file.\n\n";
                return false;
            }
            if (pixel->QueryIntAttribute("y", &y) != tinyxml2::XML_SUCCESS) {
                cerr << "Missing or invalid \'y\' attribute in config file.\n\n";
                return false;
            }
            if (fExt2BadPix.count(hdu)==0)  {
                fExt2BadPix[hdu] = std::set<std::pair<int,int>>();
            }
            fExt2BadPix[hdu].insert(make_pair(x,y));
            pixel = pixel->NextSiblingElement("pixel");
        }
    }

    /* badCols */
    tinyxml2::XMLElement* badCols = config->FirstChildElement("badCols");
    if( badCols == NULL ){
        cerr << "Missing \'badCols\' in config file.\n\n";
    } else {
        tinyxml2::XMLElement* column = badCols->FirstChildElement();
        while (column != NULL) {
            int hdu;
            if (column->QueryIntAttribute("hdu", &hdu) != tinyxml2::XML_SUCCESS) {
                cerr << "Missing or invalid \'hdu\' attribute in config file.\n\n";
                return false;
            }
            if (fExt2BadCol.count(hdu)==0)  {
                fExt2BadCol[hdu] = std::set<int>();
            }
            if (strcmp(column->Name(), "column")==0) {
                int x;
                if (column->QueryIntAttribute("x", &x) != tinyxml2::XML_SUCCESS) {
                    cerr << "Missing or invalid \'x\' attribute in config file.\n\n";
                    return false;
                }
                fExt2BadCol[hdu].insert(x);
            } else if (strcmp(column->Name(), "cRange")==0) {
                int x1, x2;
                if (column->QueryIntAttribute("x1", &x1) != tinyxml2::XML_SUCCESS) {
                    cerr << "Missing or invalid \'x1\' attribute in config file.\n\n";
                    return false;
                }
                if (column->QueryIntAttribute("x2", &x2) != tinyxml2::XML_SUCCESS) {
                    cerr << "Missing or invalid \'x2\' attribute in config file.\n\n";
                    return false;
                }
                for (int x=x1; x<=x2; x++) {
                    fExt2BadCol[hdu].insert(x);
                }
            } else {
                cerr << "Invalid element " << column->Name() << "in \'badCols\' section of config file.\n\n";
            }
            column = column->NextSiblingElement();
        }
    }

    /* bleedCols */
    tinyxml2::XMLElement* bleedCols = config->FirstChildElement("bleedCols");
    if( bleedCols == NULL ){
        cerr << "Missing \'bleedCols\' in config file.\n\n";
    } else {
        tinyxml2::XMLElement* column = bleedCols->FirstChildElement("column");
        while (column != NULL) {
            int hdu,x;
            if (column->QueryIntAttribute("hdu", &hdu) != tinyxml2::XML_SUCCESS) {
                cerr << "Missing or invalid \'hdu\' attribute in config file.\n\n";
                return false;
            }
            if (column->QueryIntAttribute("x", &x) != tinyxml2::XML_SUCCESS) {
                cerr << "Missing or invalid \'x\' attribute in config file.\n\n";
                return false;
            }
            if (fExt2BleedCol.count(hdu)==0)  {
                fExt2BleedCol[hdu] = std::set<int>();
            }
            fExt2BleedCol[hdu].insert(x);
            column = column->NextSiblingElement("column");
        }
    }

    /* bleedXEdge */
    tinyxml2::XMLElement* bleedXEdges = config->FirstChildElement("bleedXEdges");
    if( bleedXEdges == NULL ){
        cerr << "Missing \'bleedXEdges\' in config file.\n\n";
    } else {
        tinyxml2::XMLElement* edge = bleedXEdges->FirstChildElement("edge");
        while (edge != NULL) {
            int hdu,x;
            if (edge->QueryIntAttribute("hdu", &hdu) != tinyxml2::XML_SUCCESS) {
                cerr << "Missing or invalid \'hdu\' attribute in config file.\n\n";
                return false;
            }
            if (edge->QueryIntAttribute("x", &x) != tinyxml2::XML_SUCCESS) {
                cerr << "Missing or invalid \'x\' attribute in config file.\n\n";
                return false;
            }
            if (fExt2BleedXEdge.count(hdu)==0)  {
                fExt2BleedXEdge[hdu] = x;
            } else {
                cerr << "Multiple \'edge\' elements defined for a single HDU.\n\n";
                return false;
            }
            edge = edge->NextSiblingElement("edge");
        }
    }


    /* Cuts */
    if( config->FirstChildElement("cuts") == 0 ){
        cerr << "Missing \'cuts\' in config file.\n\n";
        return false;
    }

    //these are optional settings; if not defined in the XML, the default value in the constructor will be used
    config->FirstChildElement("cuts")->QueryFloatAttribute("epix", &fEpixCut);
    config->FirstChildElement("cuts")->QueryIntAttribute("bleedX", &fBleedXCut);
    config->FirstChildElement("cuts")->QueryIntAttribute("bleedY", &fBleedYCut);
    config->FirstChildElement("cuts")->QueryIntAttribute("halo", &fHaloCut);
    config->FirstChildElement("cuts")->QueryIntAttribute("edgeCut", &fEdgeCut);
    config->FirstChildElement("cuts")->QueryIntAttribute("clusterThr", &fClusterThr);
    config->FirstChildElement("cuts")->QueryIntAttribute("clusterCut", &fClusterCut);
    config->FirstChildElement("cuts")->QueryIntAttribute("clusterIgnore", &fClusterIgnore);
    config->FirstChildElement("cuts")->QueryIntAttribute("crosstalk", &fCrosstalkThr);

    config->FirstChildElement("extra")->QueryBoolAttribute("looseCluster", &fLooseCluster);
    config->FirstChildElement("extra")->QueryBoolAttribute("looseClustersNoNeighbors", &fLooseClustersNoNeighbors);
    config->FirstChildElement("extra")->QueryBoolAttribute("useDiagonalPix", &fUseDiagonalPix);
    config->FirstChildElement("extra")->QueryBoolAttribute("saveTracks", &fSaveTracks);
    config->FirstChildElement("extra")->QueryBoolAttribute("maskSmall", &fMaskSmall);

    /* Extra */
    if( config->FirstChildElement("extra") == 0 ){
        cerr << "Missing \'extra\' in config file.\n\n";
        return false;
    }

    fTracksCuts = config->FirstChildElement("extra")->Attribute("saveTrackCuts");

    /* Variables that will be saved in the NTuple */
    for(int n=0;n<gNBaseTNtupleVars;++n){
        if(n>0) fNTupleVars += ":";
        fNTupleVars += gBaseTNtupleVars[n];
    }

    fIsValid = true;
    return true;
}
