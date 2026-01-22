#ifndef ANOMALOUS_JUNCTION_COUNT_H
#define ANOMALOUS_JUNCTION_COUNT_H
//----------------------------------------------------------------
#include <map>
#include <set>
#include "bam_file.h"
#include "gtf_file.h"
#include "interval_tree.h"
//----------------------------------------------------------------
using namespace std;
//----------------------------------------------------------------
class IntervalTreeEntry
{
private:
public:

    bool isLast;
    int32_t featureIndex;
    GTFTranscript * gtfTranscript;

    IntervalTreeEntry(bool isLast,int32_t featureIndex,GTFTranscript * gtfTranscript) : isLast(isLast),featureIndex(featureIndex),gtfTranscript(gtfTranscript){}
};
//----------------------------------------------------------------
class Segment
{
private:
public:

    bool isFirst;
    bool isLast;
    bool startIsSplit;
    bool stopIsSplit;
    int32_t start;
    int32_t stop;

    Segment(bool isFirst,bool isLast,bool startIsSplit,bool stopIsSplit,int32_t start,int32_t stop) : isFirst(isFirst),isLast(isLast),startIsSplit(startIsSplit),stopIsSplit(stopIsSplit),start(start),stop(stop){}
};
//----------------------------------------------------------------
enum SegmentStatus {INTERGENIC=0,INTRONIC=1,MULTI_OVERLAP=2,START_IS_ALT=3,STOP_IS_ALT=4,BOTH_ARE_ALT=5,MATCH=6};
//----------------------------------------------------------------
class OverlapEntry
{
private:
public:

    SegmentStatus segmentStatus;
    int32_t featureIndex;

    OverlapEntry(void) : segmentStatus(INTERGENIC),featureIndex(-2){}   //-2==no_feature -1==intronic 0...==exon#
    inline void Set(SegmentStatus segmentStatus,int32_t featureIndex) { this->segmentStatus=segmentStatus; this->featureIndex=featureIndex;}
};
//----------------------------------------------------------------
class SpanningJunction
{
private:
public:

    uint32_t count;
    int32_t start1;
    int32_t stop1;
    int32_t start2;
    int32_t stop2;

    SpanningJunction(int32_t start1,int32_t stop1,int32_t start2,int32_t stop2) : count(1),start1(start1),stop1(stop1),start2(start2),stop2(stop2){}

    inline bool Merge(const SpanningJunction & junction)
    {
        if(junction.start1>stop1 || junction.stop1<start1 || junction.start2>stop2 || junction.stop2<start2) return false;

        count+=junction.count;
        start1=min(start1,junction.start1);
        stop1=max(stop1,junction.stop1);
        start2=min(start2,junction.start2);
        stop2=max(stop2,junction.stop2);

        return true;
    }
};
//----------------------------------------------------------------
class SplitJunction
{
private:
public:

    uint32_t splitCount;
    uint32_t spanningCount;
    int32_t start1;
    int32_t stop2;

    SplitJunction(void){}

    inline void Init(int32_t start1,int32_t stop2)
    {
        splitCount=1;
        spanningCount=0;
        this->start1=start1;
        this->stop2=stop2;
    }

    inline void Init(const SplitJunction & junction)
    {
        splitCount=junction.splitCount;
        spanningCount=junction.spanningCount;
        start1=junction.start1;
        stop2=junction.stop2;
    }

    inline void Merge(const SplitJunction & junction)
    {
        splitCount++;
        start1=min(start1,junction.start1);
        stop2=max(stop2,junction.stop2);
    }

    inline void Merge(int32_t start1,int32_t stop2)
    {
        splitCount++;
        this->start1=min(this->start1,start1);
        this->stop2=max(this->stop2,stop2);
    }

    inline void Init(const SpanningJunction & junction)
    {
        splitCount=0;
        spanningCount=junction.count;
        start1=junction.start1;
        stop2=junction.stop2;
    }

    inline void Merge(const SpanningJunction & junction)
    {
        spanningCount++;
        //start1=min(start1,junction.start1);
        //stop2=max(stop2,junction.stop2);
    }
};
//----------------------------------------------------------------
class KnownJunction
{
private:
public:

    uint32_t count;
    set<const GTFTranscript *> transcripts;

    KnownJunction(void) : count(0) {}
};
//----------------------------------------------------------------
class BedPEEntry
{
private:
public:

    int32_t tid1,start1,stop1,tid2,start2,stop2;
    uint32_t splitCount,spanningCount;
};
//----------------------------------------------------------------
enum FragmentStatus {DROP=0,SPLIT=1,SPAN=2};
//----------------------------------------------------------------
class AnomalousJunctionCount
{
private:

    bool countDuplicates,countSupplementary;

    int32_t maxFragLen;
    int32_t minJunctionCount;
    uint8_t minMappingQuality;
    int32_t transcriptExtension;

    string inputBamFilename;
    string gtfFilename;
    string outputPrefix;

    GTFFile gtfFile;

    vector<map<int32_t,string> > geneBoundaries;
    vector<IntervalTree<int32_t,IntervalTreeEntry> > intervalTree;
    vector<map<int32_t,map<int32_t,SplitJunction> > > splitJunctions;
    map<int32_t,map<int32_t,vector<SpanningJunction> > > interchromosomalJunctions;

    void ConstructIntervalTree(const bam_hdr_t * bamHeader,vector<GTFTranscript> & transcripts);
    static bool CreateSegments(const bam1_t * read1,const bam1_t * read2,size_t & nSegments1,size_t & nSegments2,vector<Segment> & segments);
    void DetermineSegmentOverlapStatus(const Segment & segment,const Interval<int32_t,IntervalTreeEntry> & feature,OverlapEntry & overlapEntry);
    FragmentStatus CreateSplitJunctions(size_t nSegments1,size_t nSegments2,const vector<Segment> & segments,IntervalTree<int32_t,IntervalTreeEntry> & intervalTree,map<int32_t,map<int32_t,SplitJunction> > & splitJunctions);
    static void MergeInterchromosomalJunctions(SpanningJunction & newJunction,vector<SpanningJunction> & junctions);
    void ProcessSuplementary(const bam_hdr_t * header,const bam1_t * read,BamFile & evidenceBamFile);
    bool FirstPass(void);
    static bool AddSpanToSplitJunctions(int32_t remainder,const SpanningJunction & spanningJunction,map<int32_t,map<int32_t,SplitJunction> > & splitJunctions);
    static void MergeSpanningJunctions(SpanningJunction & newJunction,multimap<int32_t,SpanningJunction> & junctions);
    void WriteBedPEEntry(const bam_hdr_t * bamHeader,ofstream & bedPEFile,BedPEEntry & bedPEEntry);
    bool WriteBedPEFile(const bam_hdr_t * bamHeader);
    bool SecondPass(void);

public:

    AnomalousJunctionCount(void);
    int Run(int argc,char * argv[]);
};
//----------------------------------------------------------------
#endif // ALGORITHM_H

