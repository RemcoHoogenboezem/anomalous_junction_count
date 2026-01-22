//----------------------------------------------------------------
// Name        : gtf_file.cpp
// Author      : Remco Hoogenboezem
// Version     :
// Copyright   :
// Description :
//----------------------------------------------------------------
#include <fstream>
#include <iostream>
#include <cstring>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include "gtf_file.h"
//----------------------------------------------------------------
class TMPTranscript //Temporary GTF transcript to sort and merge
{
private:
public:

    char strand;
    string geneName;
    string seqName;
    map<int32_t,int32_t> features;

    TMPTranscript(void){}

    void Init(const char * geneName,const char * seqName,int32_t start,int32_t stop,char strand)
    {
        this->strand=strand;
        this->geneName.assign(geneName);
        this->seqName.assign(seqName);
        features[start]=stop;
    }

    bool Add(const char * geneName,const char * seqName,int32_t start,int32_t stop,char strand)
    {
        if(this->geneName.compare(geneName)!=0)
        {
            cerr << "Error: Gene names per transcript must agree" << endl;
            return false;
        }

        if(this->seqName.compare(seqName)!=0)
        {
            cerr << "Error: Sequence names per transcript must agree" << endl;
            return false;
        }

        if(this->strand!=strand)
        {
            cerr << "Error: Strand directions per transcript must agree" << endl;
            return false;
        }

        features[start]=max(features[start],stop);   //If this exact begin position already exists keep the one with the highest end coordinate
        return true;
    }

    void MergeAndCopy(const string & transcriptID,GTFTranscript & transcript)
    {
        transcript.transcriptID=transcriptID;
        transcript.geneName=geneName;
        transcript.seqName=seqName;
        transcript.strand=strand;

        auto feature=features.begin();
        auto featureEnd=features.end();
        int32_t start=feature->first;
        int32_t stop=feature->second;

        transcript.features.clear();

        for(feature++;feature!=featureEnd;feature++)
        {
            if(feature->first<=stop+1)
            {
                stop=max(stop,feature->second);
                continue;
            }

            transcript.features.push_back(Feature(start,stop));
            start=feature->first;
            stop=feature->second;
        }

        transcript.features.push_back(Feature(start,stop));
    }
};
//----------------------------------------------------------------
char * GTFFile::trim(char * str,const char * set)
{
    char *src,*dst;

    for(src=dst=str ; *src!='\0' ; src++)
    {
        *dst=*src;

        for(const char * c=set ; *c!='\0' ; c++)
        {
            if(*dst==*c) goto NEXT;
        }

        dst++;

    NEXT:;
    }

    *dst='\0';

    return str;
}
//----------------------------------------------------------------
bool GTFFile::Open(const string & filename)
{
    //----------------------------------------------------------------
    //Open GTF file
    //----------------------------------------------------------------

    FILE * gtfFile=popen((string("bzcat -f ")+filename+string(" | zcat -f")).c_str(),"r");

    if(gtfFile==nullptr)
    {
        cerr << "Error: Could not open gtf file: " << filename << endl;
        return false;
    }

    //----------------------------------------------------------------
    //Read the file
    //----------------------------------------------------------------

    #define MAX_LINE 16384
    char line[MAX_LINE];

    unordered_map<string,TMPTranscript> transcripts;

    for(bool good=fgets(line,MAX_LINE,gtfFile)!=nullptr ; good==true ; good=fgets(line,MAX_LINE,gtfFile))
    {
        //Skip empty and comment lines

        if(line[0]=='\0' || line[0]=='#') continue;

        //Tokenize line

        char * pLine=line;
        char * pSeqName=strsep(&pLine,"\t"); strsep(&pLine,"\t"); char * pFeature=strsep(&pLine,"\t"); char * pStart=strsep(&pLine,"\t"); char * pStop=strsep(&pLine,"\t"); strsep(&pLine,"\t"); char * pStrand=strsep(&pLine,"\t"); strsep(&pLine,"\t"); char * pGroup=strsep(&pLine,"\t");

        if(pLine!=nullptr || pGroup==nullptr)
        {
            cerr << "Error: Invalid number of fields in GTF file" << endl;
            pclose(gtfFile);
            return false;
        }

        //Construct transcripts only include exon and UTR blocks

        if(strcasestr(pFeature,"exon")!=nullptr || strcasestr(pFeature,"UTR")!=nullptr)
        {
            int32_t start=atoi(pStart)-1;  //Coordinates are zero based with closed intervals
            int32_t stop=atoi(pStop)-1;

            if(stop<start)
            {
                cerr << "Error: Invalid interval detected end position can't be smaller compared to start position" << endl;
                pclose(gtfFile);
                return false;
            }

            char strand=pStrand[0];

            if((strand!='+' && strand!='-') || pStrand[1]!='\0')
            {
                cerr << "Warning: Strand field found in exon/UTR entry should be either + or -" << endl;
                continue;
            }

            char * pTranscriptID=strcasestr(pGroup,"transcript_id");

            if(pTranscriptID==nullptr)
            {
                cerr << "Warning: No \"transcript_id\" found in exon/UTR entry" << endl;
                continue;
            }

            pTranscriptID+=sizeof("transcript_id")-1;

            char * pGeneName=strcasestr(pGroup,"gene_name ");

            if(pGeneName==nullptr)
            {
                pGeneName=strcasestr(pGroup,"gene ");

                if(pGeneName==nullptr)
                {
                    cerr << "Warning: No \"gene_name\" or \"gene\" found in exon/UTR entry" << endl;
                    continue;
                }

                pGeneName+=sizeof("gene")-1;
            }
            else
            {
                pGeneName+=sizeof("gene_name")-1;
            }

            pTranscriptID=trim(strsep(&pTranscriptID,";")," \"");
            pGeneName=trim(strsep(&pGeneName,";")," \"");

            if(transcripts.count(pTranscriptID)==0)
            {
                transcripts[pTranscriptID].Init(pGeneName,pSeqName,start,stop,strand);
                continue;
            }

            if(transcripts[pTranscriptID].Add(pGeneName,pSeqName,start,stop,strand)==false)
            {
                pclose(gtfFile);
                return false;
            }
        }
    }

    pclose(gtfFile);

    if(transcripts.empty())
    {
        cerr << "Error: Could not open/read gtf file or the file wass empty: " << filename << endl;
        return false;
    }

    //----------------------------------------------------------------
    //Merge transcripts (Only keep unique transcripts)
    //----------------------------------------------------------------

    this->transcripts.clear();

    {
        GTFTranscript gtfTranscript;
        unordered_set<GTFTranscript> exists;

        for(auto & transcript : transcripts)
        {
            transcript.second.MergeAndCopy(transcript.first,gtfTranscript);

            if(exists.count(gtfTranscript)==0)
            {
                exists.insert(gtfTranscript);
                this->transcripts.push_back(gtfTranscript);
            }
        }
    }

    //----------------------------------------------------------------
    //Done
    //----------------------------------------------------------------

    return true;
}
//----------------------------------------------------------------

/*
    //----------------------------------------------------------------
    //Write bed file to test
    //----------------------------------------------------------------

    ofstream bedFile("/home/remco/Desktop/transcripts_combined.bed");

    for(const auto & transcript : this->transcripts)
    {
        int32_t chromStart=transcript.features.begin()->start;
        int32_t chromEnd=transcript.features.rbegin()->stop;
        size_t blockCount=transcript.features.size();

        bedFile << transcript.seqName << '\t' << chromStart << '\t' << chromEnd << '\t' << transcript.transcriptID << "\t0\t" << transcript.strand << '\t' << chromStart << '\t' << chromEnd << "\t255,0,0\t" << blockCount;

        bedFile << '\t' << 1+transcript.features[0].stop-transcript.features[0].start;

        for(size_t i=1;i<blockCount;i++)
        {
            bedFile << ',' << 1+transcript.features[i].stop-transcript.features[i].start;
        }

        bedFile << '\t' << transcript.features[0].start-chromStart;

        for(size_t i=1;i<blockCount;i++)
        {
            bedFile << ',' << transcript.features[i].start-chromStart;
        }

        bedFile << '\n';
    }

    bedFile.close();
*/
