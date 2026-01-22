//----------------------------------------------------------------
#ifndef BAM_FILE_H
#define	BAM_FILE_H
//----------------------------------------------------------------
#include <string>
#include "sam.h"
//----------------------------------------------------------------
using namespace std;
//----------------------------------------------------------------
class PileupData
{
private:
public:

    int nAlignments,targetId,pos;
    const bam_pileup1_t * alignments;
};
//----------------------------------------------------------------
class PileupData64
{
private:
public:

    int nAlignments,targetId;
    hts_pos_t pos;
    const bam_pileup1_t * alignments;
};
//----------------------------------------------------------------
class BamFile
{
private:

    bool canRead;
    bool canWrite;

    uint16_t requiredFlags;
    uint16_t filterFlags;

    int maxPileupDepth;

    uint64_t nTotalAlignments;

    htsFile * bamFile;
    bam_hdr_t * bamHeader;
    bam1_t * alignment;
    hts_idx_t * bamIndex;
    hts_itr_t * bamIterator;
    bam_plp_t pileupIterator;

    static int Read(void * data,bam1_t * alignment);
    static int ReadRegion(void *data, bam1_t * alignment);

public:

    BamFile(void);
    ~BamFile(void);

    const bam_hdr_t * GetHeader(void);
    const bam1_t * GetAlignment(void);
    uint64_t GetTotalAlignments(void);

    bool Open(const string & filename,uint16_t requiredFlags=0,uint16_t filterFlags=0,int maxDepth=10000);
    bool Create(const string & filename,const bam_hdr_t * bamHeader);
    void Close(void);

    inline int ReadNoCheck(void) {return sam_read1(bamFile,bamHeader,alignment);}
    inline int WriteNoCheck(const bam1_t * alignment) {return sam_write1(bamFile,bamHeader,alignment);}

    int Read(void);
    int Write(const bam1_t * alignment);

    bool SetRegion(int tid,hts_pos_t beg,hts_pos_t end);
    bool SetRegion(const string & region);
    int ReadRegion(void);

    int Pileup(PileupData & pileupData);
    int Pileup64(PileupData64 & pileupData);

    int PileupRegion(PileupData & pileupData);
    int PileupRegion64(PileupData64 & pileupData);
};
//----------------------------------------------------------------
#endif // BAM_FILE_H
