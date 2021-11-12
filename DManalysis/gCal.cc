#include "tinyxml2.h"
#include "gCal.h"
#include "globalConstants.h"
#include <iostream>
#include <sstream>


using namespace std;

gCal::gCal(): fDefaultCal(-1), fIsValid(false)
{}

float gCal::getExtCal(const int ext){

    if (fExt2Cal.count(ext)==0) {
        return fDefaultCal;
    }
    return fExt2Cal[ext];
}

bool gCal::isCrossTalk(const int ohdu, const int x, const int y){
    if (fExt2CrossTalk[1].count(make_pair(x,y))!=0 || fExt2CrossTalk[2].count(make_pair(x,y))!=0 || fExt2CrossTalk[3].count(make_pair(x,y))!=0 || fExt2CrossTalk[4].count(make_pair(x,y))!=0)
    {
        return true;
    } else {
        return false;
    }
}


void gCal::printVariables(){
    cout << endl;
    cout << "===== Configuration parameters =====\n";

    cout << "cal:\n";
    cout << "\tdefault: " << fDefaultCal << endl;
    typedef std::map<int, float>::iterator it_type;
    for(it_type it = fExt2Cal.begin(); it != fExt2Cal.end(); it++) {
        cout << "\text" << it->first << ": " << it->second << endl;
    }

    cout << "crossTalk:\n";
    typedef std::map<int, std::set<std::pair<int,int>>>::iterator it_crosstalk;
    for(it_crosstalk it = fExt2CrossTalk.begin(); it != fExt2CrossTalk.end(); it++) {
        typedef std::set<std::pair<int,int>>::iterator it_type;
        for(it_type it2 = it->second.begin(); it2 != it->second.end(); it2++) {
            cout << "\text" << it->first << " x: " << it2->first << " y: " << it2->second << endl;
        }
    }

    cout << "====================================\n";
}

bool gCal::readConfFile(const char* confFileName = "extractCal.xml"){
    fFilename = confFileName;

    tinyxml2::XMLDocument doc;
    if(doc.LoadFile( confFileName ) != 0){
        cerr << red << "\nCan't read config file! Will not continue.\n\n" << normal;
        return false;
    }
    //   doc.Print();

    tinyxml2::XMLElement* config = doc.FirstChildElement("extractCal");
    if ( config == NULL ){
        cerr << "Missing \'extractCal\' in config file.\n\n";
        return false;
    }

    tinyxml2::XMLElement* calibration = config->FirstChildElement("calibration");
    /* Calibration */
    if( calibration == NULL ){
        cerr << "Missing \'calibration\' in config file.\n\n";
        return false;
    }

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
    
    tinyxml2::XMLElement* crossTalk = config->FirstChildElement("crossTalk");
    if( crossTalk == NULL ){
        cerr << "No \'Intermodular crossTalk\' in config file.\n\n";
    } else {
        tinyxml2::XMLElement* pixel = crossTalk->FirstChildElement("pixel");
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
            if (fExt2CrossTalk.count(hdu)==0)  {
                fExt2CrossTalk[hdu] = std::set<std::pair<int,int>>();
            }
            fExt2CrossTalk[hdu].insert(make_pair(x,y));
            // cout << hdu << endl;
            // cout << x << endl;
            // cout << y << endl << endl;
            pixel = pixel->NextSiblingElement("pixel");
        }

    }


    fIsValid = true;

    return true;
}
