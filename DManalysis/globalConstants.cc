#include "globalConstants.h"

int  gVerbosity      = true;
bool gIntermodularCrossTalk = false;
const char *gBaseTNtupleVars[]  = {"runID", "ohdu", "nSat", "flag", "xMin", "xMax", "yMin", "yMax"};
const int   gNBaseTNtupleVars   = 8;
const char *gExtraTNtupleVars[] = {"e", "n", "xBary", "yBary", "xVar", "yVar", "eRaw","xy1Var","xy2Var","alpha","R","xy1Width","xy2Width","xCen","yCen"}; //xy1Var is in the axis where x*y>0 and xy2Var where x*y<0. Alpha is the angle of the main axis of the event (Y/X).
const int   gNExtraTNtupleVars  = 15;

/* Internally generated generated mask tags */
const int kHasNeighbour  = 1; // if the pixel has a neighbour of >0.5e- (or cut epix)
const int kNoEmptyMask   = 2; // if the pixel has >0.5e- (or cut epix) 
const int kInBleedZone   = 4; // if the pixel may be affected by CTI
const int kNearBigEvent  = 8; // if the pixel may be affected by High Energy Events (halo mask)
const int kCrossTalk     = 16; // if the pixel may be affected by cross talk
const int kImageTooNoisy = 32; // if image is too noisy
const int kNearEdge      = 64; // if the pixel is near the edge of the CCD
const int kOverscanEvent = 128; // if the pixel is in the same row as a Serial Register hit (check code)

/* SENSEI 2020 masks*/

const int kCluster       = 256; // look code
const int kBadPix        = 512; // if the pixel is in a pixel with an excess of single pixel events, as analyzed over a certain dataset
const int kBadCol        = 1024; // 
const int kLooseCluster  = 2048; // 
const int kExtendedBleed = 4096; // if the pixel is in a row with an excess of 1e- events, as analyzed over a certain dataset
const int kSmallCluster  = 8192; // look code
const int kFullWell      = 16384; // trapped charge to the right of a full-well cluster

//colors
const char cyan[] = { 0x1b, '[', '1', ';', '3', '6', 'm', 0 };
const char magenta[] = { 0x1b, '[', '1', ';', '3', '5', 'm', 0 };
const char red[] = { 0x1b, '[', '1', ';', '3', '1', 'm', 0 };
const char green[] = { 0x1b, '[', '1', ';', '3', '2', 'm', 0 };
const char yellow[] = { 0x1b, '[', '1', ';', '3', '3', 'm', 0 };
const char blue[] = "\x1b[1;34m";

const char bold[] = "\x1b[1;39m";

const char whiteOnRed[]    = "\x1b[1;41m";
const char whiteOnGreen[]  = "\x1b[1;42m";
const char whiteOnPurple[] = "\x1b[1;45m";
const char whiteOnViolet[] = "\x1b[1;44m";
const char whiteOnBrown[]  = "\x1b[1;43m";
const char whiteOnGray[]   = "\x1b[1;47m";

const char normal[] = { 0x1b, '[', '0', ';', '3', '9', 'm', 0 };
