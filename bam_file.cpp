//----------------------------------------------------------------
// Name        : bam_file.cpp
// Author      : Remco Hoogenboezem
// Version     :
// Copyright   :
// Description :
//----------------------------------------------------------------
#include <iostream>
#include "bam_file.h"
//----------------------------------------------------------------
BamFile::BamFile(void) : canRead(false),canWrite(false),nTotalAlignments(0),bamFile(nullptr),bamHeader(nullptr),bamIndex(nullptr),alignment(nullptr),bamIterator(nullptr),pileupIterator(nullptr){}
BamFile::~BamFile(void) {Close();}
//----------------------------------------------------------------
int BamFile::Read(void * data,bam1_t * alignment)
{
    for(;;)
    {
        int ret=sam_read1(((BamFile*)data)->bamFile,((BamFile*)data)->bamHeader,alignment); if(ret<0) return ret;
        uint16_t alignmentFlag=alignment->core.flag;
        if((alignmentFlag&((BamFile*)data)->requiredFlags)==((BamFile*)data)->requiredFlags && (alignmentFlag&((BamFile*)data)->filterFlags)==0) return ret;
    }
}
//----------------------------------------------------------------
int BamFile::ReadRegion(void *data, bam1_t * alignment)
{
    for(;;)
    {
        int ret=sam_itr_next(((BamFile*)data)->bamFile,((BamFile*)data)->bamIterator,alignment); if(ret<0) return ret;
        uint16_t alignmentFlag=alignment->core.flag;
        if((alignmentFlag&((BamFile*)data)->requiredFlags)==((BamFile*)data)->requiredFlags && (alignmentFlag&((BamFile*)data)->filterFlags)==0) return ret;
    }
}
//----------------------------------------------------------------
const bam_hdr_t * BamFile::GetHeader(void) {return bamHeader;}
//----------------------------------------------------------------
const bam1_t * BamFile::GetAlignment(void) {return alignment;}
//----------------------------------------------------------------
uint64_t BamFile::GetTotalAlignments(void) {return nTotalAlignments;}
//----------------------------------------------------------------
bool BamFile::Open(const string & filename,uint16_t requiredFlags,uint16_t filterFlags,int maxPileupDepth)
{
    //----------------------------------------------------------------
    //Open bam file
    //----------------------------------------------------------------

    htsFile * bamFile=sam_open(filename.c_str(),"r");

    if(bamFile==nullptr)
    {
        cerr << "Error: Could not open bam file: " << filename << endl;
        return false;
    }

    //----------------------------------------------------------------
    //Read header
    //----------------------------------------------------------------

    bam_hdr_t * bamHeader=sam_hdr_read(bamFile);

    if(bamHeader==nullptr)
    {
        cerr << "Error: Could not read bam header: " << filename << endl;
        sam_close(bamFile);
        return false;
    }

    //----------------------------------------------------------------
    //Load index
    //----------------------------------------------------------------

    hts_idx_t * bamIndex=sam_index_load(bamFile,filename.c_str());

    if(bamIndex==nullptr)
    {
        cerr << "Error: Could not read bam index: " << filename << endl;
        sam_hdr_destroy(bamHeader);
        sam_close(bamFile);
        return false;
    }

    //----------------------------------------------------------------
    //Calculate total read count
    //----------------------------------------------------------------

    int nTargets=sam_hdr_nref(bamHeader);

    if(nTargets==-1)
    {
        cerr << "Error: Could not calculate total read count: " << filename << endl;
        hts_idx_destroy(bamIndex);
        sam_hdr_destroy(bamHeader);
        sam_close(bamFile);
        return false;
    }

    uint64_t nTotalAlignments=0;

    for(int i=0;i<nTargets;++i)
    {
        uint64_t u, v;
        hts_idx_get_stat(bamIndex,i,&u,&v);
        nTotalAlignments+=u+v;
    }

    nTotalAlignments+=hts_idx_get_n_no_coor(bamIndex);

    //----------------------------------------------------------------
    //Update object
    //----------------------------------------------------------------

    Close();

    canRead=true;
    canWrite=false;

    this->requiredFlags=requiredFlags;
    this->filterFlags=filterFlags;

    this->maxPileupDepth=maxPileupDepth;

    this->nTotalAlignments=nTotalAlignments;

    this->bamFile=bamFile;
    this->bamHeader=bamHeader;
    this->bamIndex=bamIndex;
    this->alignment=bam_init1();


    //----------------------------------------------------------------
    //Done
    //----------------------------------------------------------------

    return true;
}
//----------------------------------------------------------------
bool BamFile::Create(const string & filename,const bam_hdr_t * bamHeader)
{
    //----------------------------------------------------------------
    //Create bam file
    //----------------------------------------------------------------

    htsFile * bamFile=sam_open(filename.c_str(),"wb");

    if(bamFile==nullptr)
    {
        cerr << "Error: Could not create bam file: " << filename << endl;
        return false;
    }

    //----------------------------------------------------------------
    //Write header
    //----------------------------------------------------------------

    if(sam_hdr_write(bamFile,bamHeader)==-1)
    {
        cerr << "Error: Could not write bam header: " << filename << endl;
        sam_close(bamFile);
        return false;
    }

    //----------------------------------------------------------------
    //Update object
    //----------------------------------------------------------------

    Close();

    canRead=false;
    canWrite=true;

    this->bamFile=bamFile;
    this->bamHeader=sam_hdr_dup(bamHeader);

    //----------------------------------------------------------------
    //Done
    //----------------------------------------------------------------

    return true;
}
//----------------------------------------------------------------
void BamFile::Close(void)
{
    if(pileupIterator!=nullptr)
    {
        bam_plp_destroy(pileupIterator);
        pileupIterator=nullptr;
    }

    if(bamIterator!=nullptr)
    {
        sam_itr_destroy(bamIterator);
        bamIterator=nullptr;
    }

    if(alignment!=nullptr)
    {
        bam_destroy1(alignment);
        alignment=nullptr;
    }

    nTotalAlignments=0;

    if(bamIndex!=nullptr)
    {
        hts_idx_destroy(bamIndex);
        bamIndex=nullptr;
    }

    if(bamHeader!=nullptr)
    {
        sam_hdr_destroy(bamHeader);
        bamHeader=nullptr;
    }

    if(bamFile!=nullptr)
    {
        sam_close(bamFile);
        bamFile=nullptr;
    }
}
//----------------------------------------------------------------
int BamFile::Read(void)
{
    if(canRead==false) return -2;

    for(;;)
    {
        int ret=sam_read1(bamFile,bamHeader,alignment); if(ret<0) return ret;
        uint16_t alignmentFlag=alignment->core.flag;
        if((alignmentFlag&requiredFlags)==requiredFlags && (alignmentFlag&filterFlags)==0) return ret;
    }
}
//----------------------------------------------------------------
int BamFile::Write(const bam1_t * alignment)
{
    if(canWrite==false) return -2;
    return sam_write1(bamFile,bamHeader,alignment);
}
//----------------------------------------------------------------
bool BamFile::SetRegion(int tid,hts_pos_t beg,hts_pos_t end)
{
    if(canRead==false) return false;

    if(pileupIterator!=nullptr)
    {
        bam_plp_destroy(pileupIterator);
        pileupIterator=nullptr;
    }

    if(bamIterator!=nullptr)
    {
        sam_itr_destroy(bamIterator);
        bamIterator=nullptr;
    }

    bamIterator=sam_itr_queryi(bamIndex,tid,beg,end);

    return bamIterator!=nullptr;
}
//----------------------------------------------------------------
bool BamFile::SetRegion(const string & region)
{
    if(canRead==false) return false;

    if(pileupIterator!=nullptr)
    {
        bam_plp_destroy(pileupIterator);
        pileupIterator=nullptr;
    }

    if(bamIterator!=nullptr)
    {
        sam_itr_destroy(bamIterator);
        bamIterator=nullptr;
    }

    bamIterator=sam_itr_querys(bamIndex,bamHeader,region.c_str());

    return bamIterator!=nullptr;
}
//----------------------------------------------------------------
int BamFile::ReadRegion(void)
{
    if(canRead==false) return -2;

    for(;;)
    {
        int ret=sam_itr_next(bamFile,bamIterator,alignment); if(ret<0) return ret;
        uint16_t alignmentFlag=alignment->core.flag;
        if((alignmentFlag&requiredFlags)==requiredFlags && (alignmentFlag&filterFlags)==0) return ret;
    }
}
//----------------------------------------------------------------
int BamFile::Pileup(PileupData & pileupData)
{
    if(canRead==false) return -2;

    if(pileupIterator==nullptr)
    {
        pileupIterator=bam_plp_init(&Read,this);

        if(pileupIterator==nullptr)
        {
            return -2;
        }

        bam_plp_set_maxcnt(pileupIterator,maxPileupDepth);
    }

    pileupData.alignments=bam_plp_auto(pileupIterator,&pileupData.targetId,&pileupData.pos,&pileupData.nAlignments);

    if(pileupData.alignments==nullptr)
    {
        bam_plp_destroy(pileupIterator);
        pileupIterator=nullptr;
        return -1;
    }

    return 0;
}
//----------------------------------------------------------------
int BamFile::Pileup64(PileupData64 & pileupData)
{
    if(canRead==false) return -2;

    if(pileupIterator==nullptr)
    {
        pileupIterator=bam_plp_init(&Read,this);

        if(pileupIterator==nullptr)
        {
            return -2;
        }

        bam_plp_set_maxcnt(pileupIterator,maxPileupDepth);
    }

    pileupData.alignments=bam_plp64_auto(pileupIterator,&pileupData.targetId,&pileupData.pos,&pileupData.nAlignments);

    if(pileupData.alignments==nullptr)
    {
        bam_plp_destroy(pileupIterator);
        pileupIterator=nullptr;
        return -1;
    }

    return 0;
}
//----------------------------------------------------------------
int BamFile::PileupRegion(PileupData & pileupData)
{
    if(bamIterator==nullptr) return -2;

    if(pileupIterator==nullptr)
    {
        pileupIterator=bam_plp_init(&ReadRegion,this);

        if(pileupIterator==nullptr)
        {
            return -2;
        }

        bam_plp_set_maxcnt(pileupIterator,maxPileupDepth);
    }

    pileupData.alignments=bam_plp_auto(pileupIterator,&pileupData.targetId,&pileupData.pos,&pileupData.nAlignments);

    if(pileupData.alignments==nullptr)
    {
        bam_plp_destroy(pileupIterator);
        pileupIterator=nullptr;
        return -1;
    }

    return 0;
}
//----------------------------------------------------------------
int BamFile::PileupRegion64(PileupData64 & pileupData)
{
    if(bamIterator==nullptr) return -2;

    if(pileupIterator==nullptr)
    {
        pileupIterator=bam_plp_init(&ReadRegion,this);

        if(pileupIterator==nullptr)
        {
            return -2;
        }

        bam_plp_set_maxcnt(pileupIterator,maxPileupDepth);
    }

    pileupData.alignments=bam_plp64_auto(pileupIterator,&pileupData.targetId,&pileupData.pos,&pileupData.nAlignments);

    if(pileupData.alignments==nullptr)
    {
        bam_plp_destroy(pileupIterator);
        pileupIterator=nullptr;
        return -1;
    }

    return 0;
}
//----------------------------------------------------------------

