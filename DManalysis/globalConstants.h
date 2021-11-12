#ifndef _globalConstants_h_
#define _globalConstants_h_


extern int  gVerbosity;
extern bool gIntermodularCrossTalk;
extern const char *gBaseTNtupleVars[];
extern const int   gNBaseTNtupleVars;
extern const char *gExtraTNtupleVars[];
extern const int   gNExtraTNtupleVars;

/* Internally generated generated mask tags */
extern const int kHasNeighbour;
extern const int kNoEmptyMask;
extern const int kInBleedZone;
extern const int kNearBigEvent;
extern const int kCrossTalk;
extern const int kImageTooNoisy;
extern const int kNearEdge;
extern const int kOverscanEvent;
extern const int kCluster;
extern const int kBadPix;
extern const int kBadCol;
extern const int kLooseCluster;
extern const int kExtendedBleed;
extern const int kSmallCluster;
extern const int kFullWell;

enum eMaskType { eNoMask, eExternalMask, ePartialMask, eComputeMask };

//colors
extern const char cyan[];
extern const char magenta[];
extern const char red[];
extern const char green[];
extern const char yellow[];
extern const char blue[];

extern const char bold[];

extern const char whiteOnRed[];
extern const char whiteOnGreen[];
extern const char whiteOnPurple[];
extern const char whiteOnViolet[];
extern const char whiteOnBrown[];
extern const char whiteOnGray[];

extern const char normal[];

#endif
