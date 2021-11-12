#include <map>
#include <set>
#include <string>

class gConfig
{
    public:
        static gConfig& getInstance()
        {
            static gConfig instance; // Guaranteed to be destroyed.
                                     // Instantiated on first use.
            return instance;
        }
        bool readConfFile(const char* confFileName);

        bool isBadPix(const int ohdu, const int x, const int y);
        bool isBadCol(const int ohdu, const int x);
        bool isBleedCol(const int ohdu, const int x);

        bool isPastBleedXEdge(const int ohdu, const int x);
              
        bool  hasCal() {return fHasCal;};
        bool  isValid(){return fIsValid;};
        float getExtCal(const int ext);
        float epixCut() {return fEpixCut;};
        int bleedXCut() {return fBleedXCut;};
        int bleedYCut() {return fBleedYCut;};
        int haloCut() {return fHaloCut;};
        int edgeCut() {return fEdgeCut;};
        int clusterThr() {return fClusterThr;};
        int clusterCut() {return fClusterCut;};
        int clusterIgnore() {return fClusterIgnore;};
        int crosstalkThr() {return fCrosstalkThr;};

        void setHaloCut(int haloCut) {fHaloCut = haloCut;};
        
        bool  getLooseCluster()    {return fLooseCluster;};
        bool  getLooseClustersNoNeighbors()    {return fLooseClustersNoNeighbors;};
        bool  getUseDiagonalPix(){return fUseDiagonalPix;};
        bool  getSaveTracks()    {return fSaveTracks;};
        bool  getMaskSmall()    {return fMaskSmall;};
        std::string getTracksCuts(){return fTracksCuts;};
        std::string getNTupleVars(){return fNTupleVars;};
        
        void printVariables();
        std::string filename(){return fFilename;};
        
    private:
        gConfig();
        
        std::map<int, std::set<std::pair<int,int>>> fExt2BadPix;
        std::map<int, std::set<int>> fExt2BadCol;
        std::map<int, std::set<int>> fExt2BleedCol;
        std::map<int, int> fExt2BleedXEdge;

        std::map<int, float> fExt2Cal;
        float fDefaultCal;

        float fEpixCut; //threshold dividing 0 and 1e- pixels (units of e-)
        int   fBleedXCut; //bleed distance, X direction (units of pix)
        int   fBleedYCut; //bleed distance, Y direction (units of pix)
        int   fHaloCut; //halo radius (units of pix)
        int   fEdgeCut; //edge distance, typically equal to halo radius (units of pix)
        int   fClusterThr; //low-E cluster threshold (units of e-)
        int   fClusterCut; //low-E cluster radius (units of pix)
        int   fClusterIgnore; //single energy bin to ignore for low-E cluster (units of e-)
        int   fCrosstalkThr; //crosstalk threshold (units of e-)

        bool fHasCal;
        bool fIsValid;
        std::string fFilename;
        
        bool fLooseCluster; //if true, run the loose cluster mask
        bool fLooseClustersNoNeighbors; //if true, only isolated pix are used to define loose clusters; if false, a 2pix2e event will form a loose cluster
        bool fUseDiagonalPix; //diagonal pix are considered neighbors for masking and clustering
        bool fMaskSmall; //mask FITS file drops mask bits >=256

        bool fSaveTracks; //save cluster hit information in ROOT tree
        std::string fTracksCuts;
        std::string fNTupleVars;
        
        // We need to declare these two. We want to make sure they
        // are unaccessable otherwise we may accidently get copies of
        // the singleton appearing.
        gConfig(gConfig const&);        // Don't Implement
        void operator=(gConfig const&); // Don't implement
};
