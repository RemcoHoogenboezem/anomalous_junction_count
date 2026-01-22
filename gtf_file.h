//----------------------------------------------------------------
#ifndef GTF_FILE_H
#define GTF_FILE_H
//----------------------------------------------------------------
#include <vector>
#include <string>
//----------------------------------------------------------------
using namespace std;
//----------------------------------------------------------------
class Feature
{
private:
public:

    int32_t start;
    int32_t stop;
    uint32_t startCount;
    uint32_t stopCount;

    Feature(int32_t start,int32_t stop) : start(start),stop(stop),startCount(0),stopCount(0){}
};
//----------------------------------------------------------------
class GTFTranscript
{
private:
public:

    string transcriptID;
    string geneName;
    string seqName;
    char strand;
    vector<Feature> features;

    GTFTranscript(void){}

    bool operator==(const GTFTranscript & other) const
    {
        size_t nFeatures=features.size();
        if(nFeatures!=other.features.size() || seqName.compare(other.seqName)!=0) return false;

        for(size_t i=0;i<nFeatures;i++)
        {
            const auto & feature=features[i];
            const auto & otherFeature=other.features[i];
            if(feature.start!=otherFeature.start || feature.stop!=otherFeature.stop) return false;
        }

        return true;
    }
};
//----------------------------------------------------------------
template<> struct std::hash<GTFTranscript>
{
    size_t operator()(const GTFTranscript & transcript) const noexcept
    {
        hash<string> stringHasher;
        hash<int32_t> intHasher;

        size_t seed=stringHasher(transcript.seqName)+0x9e3779b9;

        for(const auto & feature : transcript.features)
        {
            seed ^= intHasher(feature.start) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            seed ^= intHasher(feature.stop) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }

        return seed;
    }
};
//----------------------------------------------------------------
class GTFFile
{
private:

    static char * trim(char * str,const char * set);

public:

    vector<GTFTranscript> transcripts;
    bool Open(const string & filename);
};
//----------------------------------------------------------------
#endif // GTF_FILE_H

