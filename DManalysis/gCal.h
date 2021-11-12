#include <map>
#include <set>
#include <string>

class gCal
{
    public:
        static gCal& getInstance()
        {
            static gCal instance; // Guaranteed to be destroyed.
                                     // Instantiated on first use.
            return instance;
        }
        bool readConfFile(const char* confFileName);
        
        bool  isValid() {return fIsValid;};
        float getExtCal(const int ext);
        
        void printVariables();
        std::string filename(){return fFilename;};
        bool isCrossTalk(const int ohdu, const int x, const int y);
        
    private:
        gCal();
        
        float fDefaultCal;
        bool fIsValid;
        std::string fFilename;

        std::map<int, float> fExt2Cal;
        std::map<int, std::set<std::pair<int,int>>> fExt2CrossTalk;
        
        // We need to declare these two. We want to make sure they
        // are unaccessable otherwise we may accidently get copies of
        // the singleton appearing.
        gCal(gCal const&);        // Don't Implement
        void operator=(gCal const&); // Don't implement
};
