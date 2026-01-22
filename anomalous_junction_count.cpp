//----------------------------------------------------------------
// Name        : anomalous_junction_count.cpp
// Author      : Remco Hoogenboezem
// Version     :
// Copyright   :
// Description :
//----------------------------------------------------------------
#include <unordered_map>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <getopt.h>
#include "progress_bar.h"
#include "anomalous_junction_count.h"
//----------------------------------------------------------------
void AnomalousJunctionCount::ConstructIntervalTree(const bam_hdr_t * bamHeader,vector<GTFTranscript> & transcripts)
{
    //----------------------------------------------------------------
    //Hash chromosome names
    //----------------------------------------------------------------

    size_t nTargets=size_t(bamHeader->n_targets);
    map<string,size_t> targetIds;

    for(size_t i=0;i<nTargets;i++) targetIds.emplace(bamHeader->target_name[i],i);

    //----------------------------------------------------------------
    //Iterate over transcripts and create an interval for each feature
    //----------------------------------------------------------------

    vector<map<int32_t,set<string> > > tempGeneBoundaries(nTargets);
    vector<vector<Interval<int32_t,IntervalTreeEntry> > > intervals(nTargets);

    for(auto & transcript : transcripts)
    {
        auto it=targetIds.find(transcript.seqName);

        if(it==targetIds.end()) continue;

        tempGeneBoundaries[it->second][transcript.features.begin()->start].insert(transcript.geneName);
        tempGeneBoundaries[it->second][transcript.features.rbegin()->stop].insert(transcript.geneName);

        auto & features=transcript.features;
        size_t iLast=features.size()-1;

        intervals[it->second].push_back(Interval<int32_t,IntervalTreeEntry>(features[0].start,features[0].stop,IntervalTreeEntry(0==iLast,0,&transcript)));

        for(size_t i=1;i<=iLast;i++)
        {
            intervals[it->second].push_back(Interval<int32_t,IntervalTreeEntry>(features[i-1].stop+1,features[i].start-1,IntervalTreeEntry(false,-1,&transcript)));
            intervals[it->second].push_back(Interval<int32_t,IntervalTreeEntry>(features[i].start,features[i].stop,IntervalTreeEntry(i==iLast,int32_t(i),&transcript)));
        }
    }

    //----------------------------------------------------------------
    //Construct interval tree
    //----------------------------------------------------------------

    geneBoundaries.assign(nTargets,map<int32_t,string>());
    intervalTree.assign(nTargets,IntervalTree<int32_t,IntervalTreeEntry>());

    for(size_t i=0;i<nTargets;i++)
    {
        auto & tempBoundariesByTarget=tempGeneBoundaries[i];
        auto & boundariesByTarget=geneBoundaries[i];

        for(auto & boundary : tempBoundariesByTarget)
        {
            auto geneName=boundary.second.begin();
            auto geneNameEnd=boundary.second.end();
            string combined=*geneName;
            for(geneName++;geneName!=geneNameEnd;geneName++) combined+=','+*geneName;
            boundariesByTarget[boundary.first]=combined;
        }

        if(intervals[i].size()==0)
        {
            cerr << "Warning: No intervals in gtf file for target: " << bamHeader->target_name[i] << endl;
            continue;
        }

        intervalTree[i]=move(intervals[i]);
    }
}
//----------------------------------------------------------------
bool AnomalousJunctionCount::CreateSegments(const bam1_t * read1,const bam1_t * read2,size_t & nSegments1,size_t & nSegments2,vector<Segment> & vSegments)
{
    map<int32_t,int32_t> segments; nSegments1=0; nSegments2=0;

    //----------------------------------------------------------------
    //Create segments read1
    //----------------------------------------------------------------

    int32_t begin,end; begin=end=int32_t(read1->core.pos);
    uint32_t * cigar=bam_get_cigar(read1);
    uint32_t * cigarEnd=cigar+read1->core.n_cigar;

    for(;cigar<cigarEnd;cigar++)
    {
        int32_t opLen=int32_t(bam_cigar_oplen(*cigar));

        switch(bam_cigar_op(*cigar))
        {
            case BAM_CMATCH: case BAM_CDEL: case BAM_CEQUAL: case BAM_CDIFF:

                end+=opLen;
                continue;

            case BAM_CREF_SKIP:

                segments[begin]=end-1; nSegments1++;
                begin=end+=opLen;
                continue;
        }
    }

    segments[begin]=end-1; nSegments1++;

    //----------------------------------------------------------------
    //Create segments read2
    //----------------------------------------------------------------

    begin=end=read2->core.pos;
    cigar=bam_get_cigar(read2);
    cigarEnd=cigar+read2->core.n_cigar;

    map<int32_t,int32_t>::iterator segment;
    map<int32_t,int32_t>::iterator lastSegment=--segments.end();

    for(;cigar<cigarEnd;cigar++)
    {
        int32_t opLen=int32_t(bam_cigar_oplen(*cigar));

        switch(bam_cigar_op(*cigar))
        {
            case BAM_CMATCH: case BAM_CDEL: case BAM_CEQUAL: case BAM_CDIFF:

                end+=opLen;
                continue;

            case BAM_CREF_SKIP:

                segment=--segments.upper_bound(begin);

                if(segment==lastSegment)
                {
                    if(begin>segment->second+1)
                    {
                        segments[begin]=end-1;nSegments2++;
                        begin=end+=opLen;
                        goto NO_OVERLAP;
                    }

                    segment->second=end-1;nSegments2++;
                    begin=end+=opLen;
                    goto NO_OVERLAP;
                }

                if(begin>segment->second+1)
                {
                    if(++segment==--segments.upper_bound(end-1))
                    {
                        if(segment==lastSegment)
                        {
                            segment->second=end-1;nSegments2++;
                            begin=end+=opLen;
                            goto NO_OVERLAP;
                        }

                        if(end-1==segment->second)
                        {
                            segment++;nSegments2++;
                            begin=end+=opLen;
                            goto OVERLAP;
                        }
                    }

                    return false;
                }

                if(end-1!=segment->second) return false;

                segment++;nSegments2++;
                begin=end+=opLen;
                goto OVERLAP;
        }
    }

    segment=--segments.upper_bound(begin);

    if(segment==lastSegment)
    {
        if(begin>segment->second+1)
        {
            segments[begin]=end-1; nSegments2++;
            goto STORE_SEGMENTS;
        }

        segment->second=max(segment->second,end-1); nSegments2++;
        goto STORE_SEGMENTS;
    }

    if(begin>segment->second+1)
    {
        if(++segment==--segments.upper_bound(end-1))
        {
            if(segment==lastSegment) segment->second=max(segment->second,end-1);

            nSegments2++;
            goto STORE_SEGMENTS;
        }

        return false;
    }

    nSegments2++;
    goto STORE_SEGMENTS;

    //----------------------------------------------------------------
    //At least overlap in the next segment
    //----------------------------------------------------------------

OVERLAP:

    for(cigar++;cigar<cigarEnd;cigar++)
    {
        int32_t opLen=int32_t(bam_cigar_oplen(*cigar));

        switch(bam_cigar_op(*cigar))
        {
            case BAM_CMATCH: case BAM_CDEL: case BAM_CEQUAL: case BAM_CDIFF:

                end+=opLen;
                continue;

            case BAM_CREF_SKIP:

                if(begin!=segment->first) return false;

                if(segment==lastSegment)
                {
                    segment->second=end-1;nSegments2++;
                    begin=end+=opLen;
                    goto NO_OVERLAP;
                }

                if(end-1!=segment->second) return false;

                segment++;nSegments2++;
                begin=end+=opLen;
                continue;
        }
    }

    if(begin!=segment->first) return false;

    if(segment==lastSegment) segment->second=max(segment->second,end-1);

    nSegments2++;
    goto STORE_SEGMENTS;

    //----------------------------------------------------------------
    //No overlap in the remaining segment(s)
    //----------------------------------------------------------------

NO_OVERLAP:

    for(cigar++;cigar<cigarEnd;cigar++)
    {
        int32_t opLen=int32_t(bam_cigar_oplen(*cigar));

        switch(bam_cigar_op(*cigar))
        {
            case BAM_CMATCH: case BAM_CDEL: case BAM_CEQUAL: case BAM_CDIFF:

                end+=opLen;
                continue;

            case BAM_CREF_SKIP:

                segments[begin]=end-1; nSegments2++;
                begin=end+=opLen;
                continue;
        }
    }

    segments[begin]=end-1; nSegments2++;

STORE_SEGMENTS:

    //----------------------------------------------------------------
    //Store segments in a vector
    //----------------------------------------------------------------

    vSegments.clear();

    size_t nSegments=segments.size();
    segment=segments.begin();

    if(nSegments==nSegments1+nSegments2)    //No overlapping segments
    {
        if(nSegments1==1)
        {
            vSegments.emplace_back(true,false,false,false,segment->first,segment->second); segment++;
        }
        else
        {
            vSegments.emplace_back(true,false,false,true,segment->first,segment->second); segment++;
            for(size_t i=2;i<nSegments1;i++) {vSegments.emplace_back(false,false,true,true,segment->first,segment->second); segment++;}
            vSegments.emplace_back(false,false,true,false,segment->first,segment->second); segment++;
        }

        if(nSegments2==1)
        {
            vSegments.emplace_back(false,true,false,false,segment->first,segment->second);
        }
        else
        {
            vSegments.emplace_back(false,false,false,true,segment->first,segment->second); segment++;
            for(size_t i=2;i<nSegments2;i++) {vSegments.emplace_back(false,false,true,true,segment->first,segment->second); segment++;}
            vSegments.emplace_back(false,true,true,false,segment->first,segment->second);
        }
    }
    else //Overlapping segments
    {
        if(nSegments==1)
        {
            vSegments.emplace_back(true,true,false,false,segment->first,segment->second);
        }
        else
        {
            vSegments.emplace_back(true,false,false,true,segment->first,segment->second); segment++;
            for(size_t i=2;i<nSegments;i++) {vSegments.emplace_back(false,false,true,true,segment->first,segment->second); segment++;}
            vSegments.emplace_back(false,true,true,false,segment->first,segment->second);
        }
    }

    //----------------------------------------------------------------
    //Done
    //----------------------------------------------------------------

    return true;
}
//----------------------------------------------------------------
void AnomalousJunctionCount::DetermineSegmentOverlapStatus(const Segment & segment,const Interval<int32_t,IntervalTreeEntry> & feature,OverlapEntry & overlapEntry)
{
    overlapEntry.featureIndex=feature.value.featureIndex;

    if(segment.start==feature.start || segment.startIsSplit==false && (feature.value.featureIndex>0 || segment.start>=(feature.start-transcriptExtension)))
    {
        if(segment.stop==feature.stop || segment.stopIsSplit==false && (feature.value.isLast==false || segment.stop<=(feature.stop+transcriptExtension)))
        {
            overlapEntry.segmentStatus=MATCH;
            return;
        }

        overlapEntry.segmentStatus=STOP_IS_ALT;
        return;
    }

    if(segment.stop==feature.stop || segment.stopIsSplit==false && (feature.value.isLast==false || segment.stop<(feature.stop+transcriptExtension)))
    {
        overlapEntry.segmentStatus=START_IS_ALT;
        return;
    }

    overlapEntry.segmentStatus=BOTH_ARE_ALT;
    return;
}
//----------------------------------------------------------------
#define SPLIT_TOLERANCE                 1
#define ALTERNATIVE_START_STOP_PENALTY  1
#define SPANNING_JUNCTION_PENALTY       2
#define SPLIT_JUNCTION_PENALTY          3
//----------------------------------------------------------------
const bool segmentStatusIsAlt[][8]
{
    {false,true,true,true,true,true,true,0},    //Add 1 column for possible easier addressing
    {true,false,true,true,true,true,true,0},
    {true,true,true,true,true,true,true,0},
    {true,true,true,true,false,true,false,0},
    {true,true,true,true,true,true,true,0},
    {true,true,true,true,true,true,true,0},
    {true,true,true,true,false,true,false,0}
};
//----------------------------------------------------------------
FragmentStatus AnomalousJunctionCount::CreateSplitJunctions(size_t nSegments1,size_t nSegments2,const vector<Segment> & segments,IntervalTree<int32_t,IntervalTreeEntry> & intervalTree,map<int32_t,map<int32_t,SplitJunction> > & splitJunctions)
{
    //----------------------------------------------------------------
    //Find all overlapping features
    //----------------------------------------------------------------

    auto features=intervalTree.findOverlapping(segments.begin()->start,segments.rbegin()->stop);
    if(features.size()==0) return DROP; //No overlapping features -> intergenic - intergenic not spanning any features

    //----------------------------------------------------------------
    //Determine max feature index per transcript
    //----------------------------------------------------------------

    map<const GTFTranscript*,int32_t> maxIndexPerTranscript;

    for(const auto & feature : features)
    {
        int32_t & maxIndex=maxIndexPerTranscript.try_emplace(feature.value.gtfTranscript,feature.value.featureIndex).first->second;
        maxIndex=max(maxIndex,feature.value.featureIndex);
    }

    //----------------------------------------------------------------
    //Determine overlap status between features and segments (Finite state machine)
    //----------------------------------------------------------------

    size_t nSegments=segments.size();
    map<GTFTranscript*,vector<OverlapEntry> > overlap;

    for(const auto & feature : features)
    {
        if(maxIndexPerTranscript[feature.value.gtfTranscript]==-1) continue; //Max feature index == -1 -> intronic - intronic not spanning any features for this transcript
        overlap.try_emplace(feature.value.gtfTranscript,vector<OverlapEntry>(nSegments));

        int32_t featureStart=feature.start;
        int32_t featureStop=feature.stop;

        for(size_t i=0;i<nSegments;i++)
        {
            const auto & segment=segments[i];

            if(segment.start>featureStop || segment.stop<featureStart) continue;  //No overlap

            auto & overlapEntry=overlap[feature.value.gtfTranscript][i];

            switch(overlapEntry.segmentStatus)
            {
                case INTERGENIC:

                    if(feature.value.featureIndex==-1) { overlapEntry.Set(INTRONIC,-1); continue;}
                    DetermineSegmentOverlapStatus(segment,feature,overlapEntry);
                    continue;

                case INTRONIC:

                    if(feature.value.featureIndex==-1) continue;
                    DetermineSegmentOverlapStatus(segment,feature,overlapEntry);
                    continue;

                case START_IS_ALT:
                case STOP_IS_ALT:
                case BOTH_ARE_ALT:
                case MATCH:

                    if(feature.value.featureIndex==-1) continue;
                    overlapEntry.segmentStatus=MULTI_OVERLAP; overlapEntry.featureIndex=-2;
                    continue;

                case MULTI_OVERLAP:

                    continue;
            }
        }
    }

    if(overlap.size()==0) return DROP; //All transcipts are completely intronic in the same intron for all transcripts don't keep this read

    //----------------------------------------------------------------
    //In the next section we are going to define the split junctions
    //----------------------------------------------------------------

    map<const GTFTranscript*,uint32_t> penaltyByTranscript;
    map<const GTFTranscript*,map<int32_t,map<int32_t,SplitJunction> > > splitJunctionsByTranscript;
    set<const GTFTranscript*> spanningByTranscript;

    //----------------------------------------------------------------
    //Non-overlapping segments
    //----------------------------------------------------------------

    if(nSegments==nSegments1+nSegments2)
    {
        int32_t initFragmentLen=0; for(size_t i=0;i<nSegments;i++) initFragmentLen+=1 + segments[i].stop - segments[i].start;  //Calculate initial fragment length (This fragment is at least this long)

        //Iterate over transcripts

        for(auto & transcript : overlap)
        {
            //Alternative start in first segment

            SegmentStatus segmentStatus=transcript.second[0].segmentStatus;

            if(segmentStatus==START_IS_ALT || segmentStatus==BOTH_ARE_ALT)
            {
                int32_t pos = transcript.first->features[0].start;
                const auto & segment = segments[0];
                splitJunctionsByTranscript[transcript.first][pos-1][pos].Init(segment.start,segment.stop); penaltyByTranscript[transcript.first]+=ALTERNATIVE_START_STOP_PENALTY;
            }

            //Alternative stop in last segment

            segmentStatus=transcript.second[nSegments-1].segmentStatus;

            if(segmentStatus==STOP_IS_ALT || segmentStatus==BOTH_ARE_ALT)
            {
                int32_t pos = transcript.first->features.rbegin()->stop;
                const auto & segment = segments[nSegments-1];
                splitJunctionsByTranscript[transcript.first][pos][pos+1].Init(segment.start,segment.stop); penaltyByTranscript[transcript.first]+=ALTERNATIVE_START_STOP_PENALTY;
            }

            const auto & overlapEntries = transcript.second;

            //Look for split junctions in read1

            if(nSegments1>1)
            {
                for(size_t i=1;i<nSegments1;i++)
                {
                    const auto & overlapEntry1=overlapEntries[i-1];
                    const auto & overlapEntry2=overlapEntries[i];

                    if(segmentStatusIsAlt[overlapEntry1.segmentStatus][overlapEntry2.segmentStatus] || (overlapEntry1.featureIndex>=0 && overlapEntry2.featureIndex>=0 && overlapEntry2.featureIndex-overlapEntry1.featureIndex!=1))
                    {
                        const auto & segment1=segments[i-1];
                        const auto & segment2=segments[i];
                        splitJunctionsByTranscript[transcript.first][segment1.stop][segment2.start].Init(segment1.start,segment2.stop); penaltyByTranscript[transcript.first]+=SPLIT_JUNCTION_PENALTY;
                        continue;
                    }

                    if(overlapEntry1.segmentStatus<=INTRONIC && overlapEntry2.segmentStatus<=INTRONIC)  //Intergenic - intergenic or intronic - intronic cases
                    {
                        const auto & segment1=segments[i-1];
                        const auto & segment2=segments[i];

                        for(const auto & feature : transcript.first->features)
                        {
                            if(feature.start>=segment1.stop && feature.stop<=segment2.start)
                            {
                                splitJunctionsByTranscript[transcript.first][segment1.stop][segment2.start].Init(segment1.start,segment2.stop); penaltyByTranscript[transcript.first]+=SPLIT_JUNCTION_PENALTY;
                                break;
                            }
                        }
                    }
                }
            }

            //Look for spanning junction between read1 and read2

            {
                const auto & overlapEntry1=overlapEntries[nSegments1-1];
                const auto & overlapEntry2=overlapEntries[nSegments1];

                if(segmentStatusIsAlt[overlapEntry1.segmentStatus][overlapEntry2.segmentStatus])    //Alternative segment overlap status
                {
                    spanningByTranscript.insert(transcript.first); penaltyByTranscript[transcript.first]+=SPANNING_JUNCTION_PENALTY;
                }
                else if(overlapEntry1.segmentStatus<=INTRONIC && overlapEntry2.segmentStatus<=INTRONIC) //Intergenic - intergenic or intronic - intronic cases
                {
                    const auto & segment1=segments[nSegments1-1];
                    const auto & segment2=segments[nSegments1];

                    for(const auto & feature : transcript.first->features)
                    {
                        if(feature.start>=segment1.stop && feature.stop<=segment2.start)
                        {
                            spanningByTranscript.insert(transcript.first); penaltyByTranscript[transcript.first]+=SPANNING_JUNCTION_PENALTY;
                            break;
                        }
                    }
                }
                else    //Calculate fragment length and see if this is anomalous
                {
                    int32_t fragmentLen=initFragmentLen;
                    uint32_t featureIndex1=overlapEntries[nSegments1-1].featureIndex;
                    uint32_t featureIndex2=overlapEntries[nSegments1].featureIndex;

                    if(featureIndex1==featureIndex2)    //Both reads in same feature
                    {
                        fragmentLen+=segments[nSegments1].start-segments[nSegments1-1].stop-1;
                    }
                    else                                //Both reads in different features
                    {
                        const auto & features=transcript.first->features;
                        fragmentLen+=features[featureIndex1].stop - segments[nSegments1-1].stop;
                        for(uint32_t i=featureIndex1+1;i<featureIndex2;i++) fragmentLen+=1 + features[i].stop - features[i].start;
                        fragmentLen+=segments[nSegments1].start - features[featureIndex2].start;
                    }

                    if(fragmentLen>maxFragLen)
                    {
                        const auto & segment1=segments[nSegments1-1];
                        const auto & segment2=segments[nSegments1];
                        spanningByTranscript.insert(transcript.first); penaltyByTranscript[transcript.first]+=SPANNING_JUNCTION_PENALTY;
                    }
                }
            }

            //Look for split junctions in read2

            if(nSegments2>1)
            {
                for(size_t i=nSegments1+1;i<nSegments;i++)
                {
                    const auto & overlapEntry1=overlapEntries[i-1];
                    const auto & overlapEntry2=overlapEntries[i];

                    if(segmentStatusIsAlt[overlapEntry1.segmentStatus][overlapEntry2.segmentStatus] || (overlapEntry1.featureIndex>=0 && overlapEntry2.featureIndex>=0 && overlapEntry2.featureIndex-overlapEntry1.featureIndex!=1))
                    {
                        const auto & segment1=segments[i-1];
                        const auto & segment2=segments[i];
                        splitJunctionsByTranscript[transcript.first][segment1.stop][segment2.start].Init(segment1.start,segment2.stop); penaltyByTranscript[transcript.first]+=SPLIT_JUNCTION_PENALTY;
                        continue;
                    }

                    if(overlapEntry1.segmentStatus<=INTRONIC && overlapEntry2.segmentStatus<=INTRONIC)  //Intergenic - intergenic or intronic - intronic cases
                    {
                        const auto & segment1=segments[i-1];
                        const auto & segment2=segments[i];

                        for(const auto & feature : transcript.first->features)
                        {
                            if(feature.start>=segment1.stop && feature.stop<=segment2.start)
                            {
                                splitJunctionsByTranscript[transcript.first][segment1.stop][segment2.start].Init(segment1.start,segment2.stop); penaltyByTranscript[transcript.first]+=SPLIT_JUNCTION_PENALTY;
                                break;
                            }
                        }
                    }
                }
            }

            if(penaltyByTranscript[transcript.first]==0)
            {
                auto & features=transcript.first->features;
                auto & overlapEntries=transcript.second;

                for(size_t i=1;i<nSegments1;i++){features[overlapEntries[i-1].featureIndex].stopCount++; features[overlapEntries[i].featureIndex].startCount++;}
                for(size_t i=nSegments1+1;i<nSegments;i++){features[overlapEntries[i-1].featureIndex].stopCount++; features[overlapEntries[i].featureIndex].startCount++;}
                return DROP;
            }
        }
    }

    //----------------------------------------------------------------
    //Overlapping segments
    //----------------------------------------------------------------

    else
    {
        for(auto & transcript : overlap)
        {
            //Alternative start in first segment

            SegmentStatus segmentStatus=transcript.second[0].segmentStatus;

            if(segmentStatus==START_IS_ALT || segmentStatus==BOTH_ARE_ALT)
            {
                int32_t pos = transcript.first->features[0].start;
                const auto & segment = segments[0];
                splitJunctionsByTranscript[transcript.first][pos-1][pos].Init(segment.start,segment.stop); penaltyByTranscript[transcript.first]+=ALTERNATIVE_START_STOP_PENALTY;
            }

            //Alternative stop in last segment

            segmentStatus=transcript.second[nSegments-1].segmentStatus;

            if(segmentStatus==STOP_IS_ALT || segmentStatus==BOTH_ARE_ALT)
            {
                int32_t pos = transcript.first->features.rbegin()->stop;
                const auto & segment = segments[nSegments-1];
                splitJunctionsByTranscript[transcript.first][pos][pos+1].Init(segment.start,segment.stop); penaltyByTranscript[transcript.first]+=ALTERNATIVE_START_STOP_PENALTY;
            }

            const auto & overlapEntries = transcript.second;

            for(size_t i=1;i<nSegments;i++)
            {
                const auto & overlapEntry1=overlapEntries[i-1];
                const auto & overlapEntry2=overlapEntries[i];

                if(segmentStatusIsAlt[overlapEntry1.segmentStatus][overlapEntry2.segmentStatus] || (overlapEntry1.featureIndex>=0 && overlapEntry2.featureIndex>=0 && overlapEntry2.featureIndex-overlapEntry1.featureIndex!=1))
                {
                    const auto & segment1=segments[i-1];
                    const auto & segment2=segments[i];
                    splitJunctionsByTranscript[transcript.first][segment1.stop][segment2.start].Init(segment1.start,segment2.stop); penaltyByTranscript[transcript.first]+=SPLIT_JUNCTION_PENALTY;
                    continue;
                }

                if(overlapEntry1.segmentStatus<=INTRONIC && overlapEntry2.segmentStatus<=INTRONIC)  //Intergenic - intergenic or intronic - intronic cases
                {
                    const auto & segment1=segments[i-1];
                    const auto & segment2=segments[i];

                    for(const auto & feature : transcript.first->features)
                    {
                        if(feature.start>=segment1.stop && feature.stop<=segment2.start)
                        {
                            splitJunctionsByTranscript[transcript.first][segment1.stop][segment2.start].Init(segment1.start,segment2.stop); penaltyByTranscript[transcript.first]+=SPLIT_JUNCTION_PENALTY;
                            break;
                        }
                    }
                }
            }

            if(penaltyByTranscript[transcript.first]==0)
            {
                auto & features=transcript.first->features;
                auto & overlapEntries=transcript.second;

                for(size_t i=1;i<nSegments;i++) {features[overlapEntries[i-1].featureIndex].stopCount++; features[overlapEntries[i].featureIndex].startCount++;}
                return DROP;
            }
        }
    }

    //----------------------------------------------------------------
    //Sort transcripts by penalty
    //----------------------------------------------------------------

    multimap<uint32_t,const GTFTranscript*> transcriptByPenalty;
    for(const auto & it : penaltyByTranscript) transcriptByPenalty.emplace(it.second,it.first);

    //----------------------------------------------------------------
    //Merge all junctions from the transcripts with lowest penalty score
    //----------------------------------------------------------------

    map<int32_t,map<int32_t,SplitJunction> > newSplitJunctions;
    bool gotSpanning=false;

    uint32_t minPenalty=transcriptByPenalty.begin()->first;

    for(const auto & transcript : transcriptByPenalty)
    {
        if(transcript.first>minPenalty) break;

        if(splitJunctionsByTranscript.count(transcript.second)!=0)
        {
            for(const auto & pos1 : splitJunctionsByTranscript[transcript.second])
            {
                for(const auto & pos2 : pos1.second)
                {
                    newSplitJunctions[pos1.first][pos2.first].Init(pos2.second);
                }
            }
        }

        if(spanningByTranscript.count(transcript.second)!=0) gotSpanning=true;
    }

    //----------------------------------------------------------------
    //Add split junctions to final list allow tolerance in either direction
    //----------------------------------------------------------------

    map<int32_t,map<int32_t,bool> > visited;

    for(const auto & newJuncPos1 : newSplitJunctions)
    {
        auto juncPos1=splitJunctions.lower_bound(newJuncPos1.first-SPLIT_TOLERANCE);

        if(juncPos1==splitJunctions.end() || abs(newJuncPos1.first-juncPos1->first)>SPLIT_TOLERANCE)
        {
            for(const auto & newJuncPos2 : newJuncPos1.second)
            {
                splitJunctions[newJuncPos1.first][newJuncPos2.first].Init(newJuncPos2.second); visited[newJuncPos1.first][newJuncPos2.first]=true;
            }
            continue;
        }

        auto & existing=juncPos1->second;

        for(const auto & newJuncPos2 : newJuncPos1.second)
        {
            auto juncPos2=existing.lower_bound(newJuncPos2.first-SPLIT_TOLERANCE);

            if(juncPos2==existing.end() || abs(newJuncPos2.first-juncPos2->first)>SPLIT_TOLERANCE)
            {
                existing[newJuncPos2.first].Init(newJuncPos2.second); visited[juncPos1->first][newJuncPos2.first]=true;
                continue;
            }

            if(visited[juncPos1->first][juncPos2->first]==false)
            {
                juncPos2->second.Merge(newJuncPos2.second); visited[juncPos1->first][juncPos2->first]=true;
            }
        }
    }

    //----------------------------------------------------------------
    //Done
    //----------------------------------------------------------------

    return gotSpanning ? SPAN : SPLIT;
}
//----------------------------------------------------------------
void AnomalousJunctionCount::MergeInterchromosomalJunctions(SpanningJunction & newJunction,vector<SpanningJunction> & junctions)
{
    auto junctionsEnd=junctions.end();

    for(auto junction=junctions.begin();junction!=junctionsEnd;junction++)
    {
        if(newJunction.Merge(*junction))
        {
            junctions.erase(junction);
            MergeInterchromosomalJunctions(newJunction,junctions);
            return;
        }
    }

    junctions.push_back(newJunction);
}
//----------------------------------------------------------------
void AnomalousJunctionCount::ProcessSuplementary(const bam_hdr_t * header,const bam1_t * read,BamFile & evidenceBamFile)
{
    int32_t readPos=read->core.pos;

    //----------------------------------------------------------------
    //Find SA tag
    //----------------------------------------------------------------

    char * SA=(char*)bam_aux_get(read,"SA");
    if(SA++==nullptr) return;

    //----------------------------------------------------------------
    //Tokenize SA tag
    //----------------------------------------------------------------

    char newSA[1024]; char * pNewSA=newSA;
    strcpy(newSA,SA);
    char * pRName=strsep(&pNewSA,",;");
    char * pPos=strsep(&pNewSA,",;");
    strsep(&pNewSA,",;");
    char * pCIGAR=strsep(&pNewSA,",;");
    strsep(&pNewSA,",;");

    if(strsep(&pNewSA,",;")==nullptr) return;

    int32_t saPos=atoi(pPos);

    //----------------------------------------------------------------
    //Same chromosome?
    //----------------------------------------------------------------

    if(strcmp(pRName,header->target_name[read->core.tid])==0)
    {
        auto & splitJunctions=this->splitJunctions[read->core.tid];

        if(readPos<saPos)
        {
            //----------------------------------------------------------------
            //Calculate last segment for read
            //----------------------------------------------------------------

            int32_t readBegin,readEnd; readBegin=readEnd=readPos;
            uint32_t * cigar=bam_get_cigar(read);
            uint32_t * cigarEnd=cigar+read->core.n_cigar;

            for(;cigar<cigarEnd;cigar++)
            {
                int32_t opLen=bam_cigar_oplen(*cigar);

                switch(bam_cigar_op(*cigar))
                {
                    case BAM_CMATCH: case BAM_CDEL: case BAM_CEQUAL: case BAM_CDIFF:

                        readEnd+=opLen;
                        continue;

                    case BAM_CREF_SKIP:

                        readBegin=readEnd+=opLen;
                        continue;
                }
            }

            readEnd--;

            //----------------------------------------------------------------
            //Calculate first SA entry segment
            //----------------------------------------------------------------

            int32_t saBegin,saEnd; saBegin=saEnd=saPos;

            while(pCIGAR[0]!='\0')
            {
                int32_t opLen=std::strtol(pCIGAR,&pCIGAR,10);
                switch(pCIGAR++[0])
                {
                    case 'M':
                    case 'D':
                    case '=':
                    case 'X':
                        saEnd+=opLen;
                        continue;
                    case 'N':
                        break;
                }
            }

            saEnd--;

            if(readEnd>saBegin) return;

            //----------------------------------------------------------------
            //Add split junction
            //----------------------------------------------------------------

            auto juncPos1=splitJunctions.lower_bound(readEnd-SPLIT_TOLERANCE);

            if(juncPos1==splitJunctions.end() || abs(readEnd-juncPos1->first)>SPLIT_TOLERANCE)
            {
                splitJunctions[readEnd][saBegin].Init(readBegin,saEnd);
            }
            else
            {
                auto & existing=juncPos1->second;
                auto juncPos2=juncPos1->second.lower_bound(saBegin-SPLIT_TOLERANCE);

                if(juncPos2==existing.end() || abs(saBegin-juncPos2->first)>SPLIT_TOLERANCE)
                {
                    existing[saBegin].Init(readBegin,saEnd);
                }
                else
                {
                    juncPos2->second.Merge(readBegin,saEnd);
                }
            }

            //----------------------------------------------------------------
            //Write read to file
            //----------------------------------------------------------------

            if(evidenceBamFile.WriteNoCheck(read)<0)
            {
                cerr << "Error: Could not write read to evidence bam file" << endl;
                return;
            }

            return;
        }

        //----------------------------------------------------------------
        //Calculate last SA entry segment
        //----------------------------------------------------------------

        int32_t saBegin,saEnd; saBegin=saEnd=saPos;

        while(pCIGAR[0]!='\0')
        {
            int32_t opLen=std::strtol(pCIGAR,&pCIGAR,10);
            switch(pCIGAR++[0])
            {
                case 'M':
                case 'D':
                case '=':
                case 'X':
                    saEnd+=opLen;
                    continue;
                case 'N':
                    saBegin=saEnd+=opLen;
                    continue;
            }
        }

        saEnd--;

        //----------------------------------------------------------------
        //Calculate last segment for read
        //----------------------------------------------------------------

        int32_t readBegin,readEnd; readBegin=readEnd=readPos;
        uint32_t * cigar=bam_get_cigar(read);
        uint32_t * cigarEnd=cigar+read->core.n_cigar;

        for(;cigar<cigarEnd;cigar++)
        {
            int32_t opLen=bam_cigar_oplen(*cigar);

            switch(bam_cigar_op(*cigar))
            {
                case BAM_CMATCH: case BAM_CDEL: case BAM_CEQUAL: case BAM_CDIFF:

                    readEnd+=opLen;
                    continue;

                case BAM_CREF_SKIP:

                    break;
            }
        }

        readEnd--;

        if(saEnd>readBegin) return;

        //----------------------------------------------------------------
        //Add split junction
        //----------------------------------------------------------------

        auto juncPos1=splitJunctions.lower_bound(saEnd-SPLIT_TOLERANCE);

        if(juncPos1==splitJunctions.end() || abs(saEnd-juncPos1->first)>SPLIT_TOLERANCE)
        {
            splitJunctions[saEnd][readBegin].Init(saBegin,readEnd);
        }
        else
        {
            auto & existing=juncPos1->second;
            auto juncPos2=juncPos1->second.lower_bound(readBegin-SPLIT_TOLERANCE);

            if(juncPos2==existing.end() || abs(readBegin-juncPos2->first)>SPLIT_TOLERANCE)
            {
                existing[readBegin].Init(saBegin,readEnd);
            }
            else
            {
                juncPos2->second.Merge(saBegin,readEnd);
            }
        }

        //----------------------------------------------------------------
        //Write read to file
        //----------------------------------------------------------------

        if(evidenceBamFile.WriteNoCheck(read)<0)
        {
            cerr << "Error: Could not write read to evidence bam file" << endl;
            return;
        }

        return;
    }

    //----------------------------------------------------------------
    //Interchromosomal
    //----------------------------------------------------------------

    int32_t readBegin=readPos;
    int32_t readEnd=bam_endpos(read)-1;
    int32_t saBegin=saPos;
    int32_t saEnd=saPos;

    while(pCIGAR[0]!='\0')
    {
        int32_t opLen=std::strtol(pCIGAR,&pCIGAR,10);
        switch(pCIGAR++[0])
        {
            case 'M':
            case 'D':
            case 'N':
            case '=':
            case 'X':
                saEnd+=opLen;
                continue;
        }
    }

    saEnd--;

    int32_t tid=0; for(;tid<header->n_targets;tid++) if(strcmp(header->target_name[tid],pRName)==0) break;
    if(tid==header->n_targets) return;

    if(read->core.tid<tid)  //Add junction in the same order as the header file
    {
        string tid1Name(header->target_name[read->core.tid]);
        string tid2Name(pRName);

        SpanningJunction interchromosomalJunction(readBegin,readEnd,saBegin,saEnd);
        MergeInterchromosomalJunctions(interchromosomalJunction,interchromosomalJunctions[read->core.tid][tid]);
    }
    else
    {
        string tid1Name(pRName);
        string tid2Name(header->target_name[read->core.tid]);

        SpanningJunction interchromosomalJunction(saBegin,saEnd,readBegin,readEnd);
        MergeInterchromosomalJunctions(interchromosomalJunction,interchromosomalJunctions[read->core.tid][tid]);
    }

    //----------------------------------------------------------------
    //Write read to file
    //----------------------------------------------------------------

    if(evidenceBamFile.WriteNoCheck(read)<0)
    {
        cerr << "Error: Could not write read to evidence bam file" << endl;
        return;
    }
}
//----------------------------------------------------------------
bool AnomalousJunctionCount::FirstPass(void)
{
    //----------------------------------------------------------------
    //Open input bam file
    //----------------------------------------------------------------

    cerr << "Info: Open bam file" << endl;

    BamFile inputBamFile;
    if(inputBamFile.Open(inputBamFilename.c_str())==false) return false;

    const bam_hdr_t * bamHeader=inputBamFile.GetHeader();

    size_t nAlignmentsTotal=inputBamFile.GetTotalAlignments();
    cerr << "Info: Number of alignments to process: " << nAlignmentsTotal << endl;

    //----------------------------------------------------------------
    //Open gtf file (GTF file is stored in object)
    //----------------------------------------------------------------

    cerr << "Info: Open GTF file" << endl;

    if(gtfFile.Open(gtfFilename)==false) return false;

    //----------------------------------------------------------------
    //Create evidence bam file
    //----------------------------------------------------------------

    cerr << "Info: Create evidence bam file" << endl;

    BamFile evidenceBamFile;
    string evidenceBamFilename=outputPrefix+"_evidence.bam";

    if(evidenceBamFile.Create(evidenceBamFilename,bamHeader)==false) return false;

    //----------------------------------------------------------------
    //Create spanning bam file
    //----------------------------------------------------------------

    cerr << "Info: Create spanning bam file" << endl;

    BamFile spanningBamFile;
    string spanningBamFilename=outputPrefix+"_spanning.bam";

    if(spanningBamFile.Create(spanningBamFilename,bamHeader)==false) return false;

    //----------------------------------------------------------------
    //Create interval tree
    //----------------------------------------------------------------

    cerr << "Info: Create interval tree" << endl;

    ConstructIntervalTree(bamHeader,gtfFile.transcripts);

    //----------------------------------------------------------------
    //Init junctions
    //----------------------------------------------------------------

    cerr << "Info: Initialize junctions" << endl;

    splitJunctions.assign(bamHeader->n_targets,map<int32_t,map<int32_t,SplitJunction> >());
    interchromosomalJunctions.clear();

    //----------------------------------------------------------------
    //Iterate through the input bam file
    //----------------------------------------------------------------

    cerr << "Info: Process input bam file" << endl;

    uint16_t excludeFlags=BAM_FUNMAP|BAM_FMUNMAP|BAM_FSECONDARY|BAM_FQCFAIL;
    if(countDuplicates==false) excludeFlags|=BAM_FDUP;
    if(countSupplementary==false) excludeFlags|=BAM_FSUPPLEMENTARY;

    size_t nAlignmentsProcessed=0;
    size_t progress=-1;
    size_t nPairs=0;
    size_t nInconsistentPairs=0;

    const bam1_t * read2=inputBamFile.GetAlignment();

    string queryName;
    vector<Segment> segments;
    unordered_map<string,bam1_t *> mates;

    while(inputBamFile.ReadNoCheck()>=0)
    {
        //----------------------------------------------------------------
        //Show progress
        //----------------------------------------------------------------

        size_t newProgress=++nAlignmentsProcessed*100 / nAlignmentsTotal;
        if(newProgress!=progress){ cerr << progressBar[newProgress] << flush; progress=newProgress;}

        //----------------------------------------------------------------
        //Only keep properly aligned reads on target IDs in database
        //----------------------------------------------------------------

        if((read2->core.flag & BAM_FPAIRED) == 0 || (read2->core.flag & excludeFlags)!=0 || (intervalTree[read2->core.tid].empty() && intervalTree[read2->core.mtid].empty())) continue;

        //----------------------------------------------------------------
        //Process supplementary alignments (No need to group and use SA tag)
        //----------------------------------------------------------------

        if(read2->core.flag & BAM_FSUPPLEMENTARY)
        {
            ProcessSuplementary(bamHeader,read2,evidenceBamFile);
            continue;
        }

        //----------------------------------------------------------------
        //Group reads
        //----------------------------------------------------------------

        queryName=string(bam_get_qname(read2))+':'+to_string(read2->core.mpos)+':'+to_string(read2->core.pos);

        if(mates.count(queryName)==0)
        {
            queryName=string(bam_get_qname(read2))+':'+to_string(read2->core.pos)+':'+to_string(read2->core.mpos);
            mates[queryName]=bam_dup1(read2);
            continue;
        }

        bam1_t * read1=mates.at(queryName); mates.erase(queryName);

        //----------------------------------------------------------------
        //Skipp fragments that have a poor alignment score on both mates
        //----------------------------------------------------------------

        if(read1->core.qual<minMappingQuality && read2->core.qual<minMappingQuality){bam_destroy1(read1);continue;}

        nPairs++;

        //----------------------------------------------------------------
        //Reads align to the same chromosome?
        //----------------------------------------------------------------

        if(read1->core.tid==read2->core.tid)
        {
            //----------------------------------------------------------------
            //Create segments
            //----------------------------------------------------------------

            size_t nSegments1,nSegments2;
            if(CreateSegments(read1,read2,nSegments1,nSegments2,segments)==false){bam_destroy1(read1);nInconsistentPairs++;continue;}

            //----------------------------------------------------------------
            //Skip reads with only one segment
            //----------------------------------------------------------------

            size_t nSegments=segments.size();
            if(nSegments==1){bam_destroy1(read1);continue;}

            //----------------------------------------------------------------
            //Skip unspliced reads with a small insert size
            //----------------------------------------------------------------

            if(nSegments==2 && nSegments1==1 && nSegments2==1 && 1+segments[1].stop-segments[0].start<=maxFragLen) {bam_destroy1(read1);continue;}

            //----------------------------------------------------------------
            //Create junctions
            //----------------------------------------------------------------

            size_t tid=size_t(read1->core.tid);

            FragmentStatus status=CreateSplitJunctions(nSegments1,nSegments2,segments,intervalTree[tid],splitJunctions[tid]);
            if(status==DROP) continue;

            if(evidenceBamFile.WriteNoCheck(read1)<0)
            {
                cerr << "Error: Could not write read to evidence bam file" << endl;
                bam_destroy1(read1);
                for(const auto & it : mates) bam_destroy1(it.second);
                return false;
            }

            if(evidenceBamFile.WriteNoCheck(read2)<0)
            {
                cerr << "Error: Could not write read to evidence bam file" << endl;
                bam_destroy1(read1);
                for(const auto & it : mates) bam_destroy1(it.second);
                return false;
            }

            if(status==SPAN)
            {
                if(spanningBamFile.WriteNoCheck(read1)<0)
                {
                    cerr << "Error: Could not write read to spanning bam file" << endl;
                    bam_destroy1(read1);
                    for(const auto & it : mates) bam_destroy1(it.second);
                    return false;
                }

                if(spanningBamFile.WriteNoCheck(read2)<0)
                {
                    cerr << "Error: Could not write read to spanning bam file" << endl;
                    bam_destroy1(read1);
                    for(const auto & it : mates) bam_destroy1(it.second);
                    return false;
                }
            }
        }

        //----------------------------------------------------------------
        //Reads align to different chromosomes?
        //----------------------------------------------------------------

        else
        {
            string tid1Name(bamHeader->target_name[read1->core.tid]);
            string tid2Name(bamHeader->target_name[read2->core.tid]);

            SpanningJunction interchromosomalJunction(read1->core.pos,bam_endpos(read1)-1,read2->core.pos,bam_endpos(read2)-1);
            MergeInterchromosomalJunctions(interchromosomalJunction,interchromosomalJunctions[read1->core.tid][read2->core.tid]);

            if(evidenceBamFile.WriteNoCheck(read1)<0)
            {
                cerr << "Error: Could not write read to evidence bam file" << endl;
                bam_destroy1(read1);
                for(const auto & it : mates) bam_destroy1(it.second);
                return false;
            }

            if(evidenceBamFile.WriteNoCheck(read2)<0)
            {
                cerr << "Error: Could not write read to evidence bam file" << endl;
                bam_destroy1(read1);
                for(const auto & it : mates) bam_destroy1(it.second);
                return false;
            }
        }

        bam_destroy1(read1);
    }

    cerr << endl;

    //----------------------------------------------------------------
    //Clean-up
    //----------------------------------------------------------------

    cerr << "Info: Clean up" << endl;

    if(mates.size())cerr << "Info: Mate buffer is not empty, still containing: " << mates.size() << " mates";
    for(const auto & it : mates) bam_destroy1(it.second);
    inputBamFile.Close();
    evidenceBamFile.Close();
    spanningBamFile.Close();

    //----------------------------------------------------------------
    //Output number of inconsistent mates
    //----------------------------------------------------------------

    cerr << fixed;
    cerr << setprecision(2);
    cerr << "Info: Number of inconsistent mates detected: " << nInconsistentPairs << " (" << double(nInconsistentPairs*100)/double(nPairs) << "%)" << endl;

    //----------------------------------------------------------------
    //Sort and index the evidence bam file using samtools (Done)
    //----------------------------------------------------------------

    cerr << "Info: Sort and index evidence file" << endl;   //Dont sort and index the spanning file

    return system((string("samtools sort --write-index -o ")+outputPrefix+string("_evidence.bam##idx##")+outputPrefix+string("_evidence.bam.bai ")+outputPrefix+string("_evidence.bam")).c_str())==0;
}
//----------------------------------------------------------------
bool AnomalousJunctionCount::AddSpanToSplitJunctions(int32_t remainder,const SpanningJunction & spanningJunction,map<int32_t,map<int32_t,SplitJunction> > & splitJunctions)
{
    int32_t stop1=spanningJunction.stop1;
    int32_t start2=spanningJunction.start2;

    auto itSplitJuncPos1End=splitJunctions.end();

    for(auto itSplitJuncPos1=splitJunctions.lower_bound(stop1);itSplitJuncPos1!=itSplitJuncPos1End;itSplitJuncPos1++)
    {
        int32_t diffPos1=itSplitJuncPos1->first-stop1;
        if(diffPos1>remainder) break;

        auto & splitJuncPos1=itSplitJuncPos1->second;
        auto itSplitJuncPos2End=splitJuncPos1.end();

        for(auto itSplitJuncPos2=splitJuncPos1.begin();itSplitJuncPos2!=itSplitJuncPos2End;itSplitJuncPos2++)
        {
            if(itSplitJuncPos2->first>start2) continue;

            int32_t diffPos2=start2-itSplitJuncPos2->first;
            if(diffPos1+diffPos2>remainder) break;

            itSplitJuncPos2->second.Merge(spanningJunction);
            return true;
        }
    }

    return false;
}
//----------------------------------------------------------------
void AnomalousJunctionCount::MergeSpanningJunctions(SpanningJunction & newJunction,multimap<int32_t,SpanningJunction> & junctions)
{
    auto junctionsEnd=junctions.end();

    for(auto junction=junctions.lower_bound(newJunction.start1);junction!=junctionsEnd;junction++)
    {
        if(newJunction.Merge(junction->second))
        {
            junctions.erase(junction);
            MergeSpanningJunctions(newJunction,junctions);
            return;
        }
    }

    junctions.insert(make_pair(newJunction.stop2,newJunction));
}
//----------------------------------------------------------------
class StrandCountPair
{
private:
public:

    char strand;
    bool isExonic;
    uint32_t count;

    StrandCountPair(void) : isExonic(false),strand('.'),count(0){}
    inline void Set(char strand,uint32_t count){this->strand=strand;this->count=count;}
};
//----------------------------------------------------------------
void AnomalousJunctionCount::WriteBedPEEntry(const bam_hdr_t * bamHeader,ofstream & bedPEFile,BedPEEntry & bedPEEntry)
{
    //----------------------------------------------------------------
    // Determine overlapping features 1
    //----------------------------------------------------------------

    bool gotExonic1=false;
    map<string,StrandCountPair> countByGene1;

    intervalTree[bedPEEntry.tid1].visit_overlapping(bedPEEntry.stop1,bedPEEntry.stop1,[&](const Interval<int32_t,IntervalTreeEntry> & entry)
    {
        int32_t featureIndex=entry.value.featureIndex;
        GTFTranscript * gtfTranscript=entry.value.gtfTranscript;
        StrandCountPair & strandCountPair=countByGene1[gtfTranscript->geneName];

        if(strandCountPair.strand!=gtfTranscript->strand)   //Conflicint strands same gene name?
        {
            if(strandCountPair.strand=='.') {strandCountPair.strand=gtfTranscript->strand;} else {strandCountPair.strand='?';}
        }

        if(featureIndex>=0) //Exonic?
        {
            gotExonic1=strandCountPair.isExonic=true;
            if(bedPEEntry.tid1!=bedPEEntry.tid2 && gtfTranscript->strand=='-') {strandCountPair.count+=gtfTranscript->features[featureIndex].startCount;return;}
            strandCountPair.count+=gtfTranscript->features[featureIndex].stopCount; return;
        }

        auto & features=gtfTranscript->features; size_t nFeatures=features.size();  //Intronic

        for(size_t i=1;i<nFeatures;i++)
        {
            if(features[i-1].stop<bedPEEntry.stop1 && bedPEEntry.stop1<features[i].start)
            {
                strandCountPair.count+=features[i-1].stopCount; return;
            }
        }
    });

    const char intronExon[][5]={"(in)","(ex)"};

    string name1;
    string strand1;

    if(countByGene1.empty()==false)
    {
        auto entry=countByGene1.begin();
        auto entryEnd=countByGene1.end();

        for(;entry!=entryEnd;entry++)
        {
            if(gotExonic1==true && entry->second.isExonic==false) continue;
            name1=entry->first+string(intronExon[entry->second.isExonic]);
            strand1=entry->second.strand;
            break;
        }

        for(entry++;entry!=entryEnd;entry++)
        {
            if(gotExonic1==true && entry->second.isExonic==false) continue;
            name1+=string(",")+entry->first+string(intronExon[entry->second.isExonic]);
            strand1+=string(",")+entry->second.strand;
        }
    }
    else
    {
        auto & boundariesByTarget=geneBoundaries[bedPEEntry.tid1];
        auto downstreamGene=boundariesByTarget.upper_bound(bedPEEntry.stop1);

        if(downstreamGene!=boundariesByTarget.end())
        {
            int32_t disToDownstream=downstreamGene->first-bedPEEntry.stop1;

            if(downstreamGene!=boundariesByTarget.begin())
            {
               auto upstreamGene=downstreamGene;--upstreamGene;

               int32_t disToUpstream=bedPEEntry.stop1-upstreamGene->first;

                if(disToUpstream<disToDownstream)
                {
                    name1=upstreamGene->second+string("(")+to_string(disToUpstream)+string(")");
                }
                else
                {
                    name1=string("(")+to_string(disToDownstream)+string(")")+downstreamGene->second;
                }
            }
            else
            {
                name1=string("(")+to_string(disToDownstream)+string(")")+downstreamGene->second;
            }
        }
        else
        {
            if(downstreamGene!=boundariesByTarget.begin())
            {
                auto upstreamGene=downstreamGene;--upstreamGene;
                name1=upstreamGene->second+string("(")+to_string(bedPEEntry.stop1-upstreamGene->first)+string(")");
            }
        }

        strand1=string(".");
    }

    //----------------------------------------------------------------
    // Determine overlapping features 2
    //----------------------------------------------------------------

    bool gotExonic2=false;
    map<string,StrandCountPair> countByGene2;

    intervalTree[bedPEEntry.tid2].visit_overlapping(bedPEEntry.start2,bedPEEntry.start2,[&](const Interval<int32_t,IntervalTreeEntry> & entry)
    {
       int32_t featureIndex=entry.value.featureIndex;
        GTFTranscript * gtfTranscript=entry.value.gtfTranscript;
        StrandCountPair & strandCountPair=countByGene2[gtfTranscript->geneName];

        if(strandCountPair.strand!=gtfTranscript->strand)
        {
            if(strandCountPair.strand=='.') {strandCountPair.strand=gtfTranscript->strand;} else {strandCountPair.strand='?';}
        }

        if(featureIndex>=0)
        {
            gotExonic2=strandCountPair.isExonic=true;
            if(bedPEEntry.tid1!=bedPEEntry.tid2 && gtfTranscript->strand=='+') {strandCountPair.count+=gtfTranscript->features[featureIndex].stopCount;return;}
            strandCountPair.count+=gtfTranscript->features[featureIndex].startCount; return;
        }

        auto & features=gtfTranscript->features; size_t nFeatures=features.size();

        for(size_t i=1;i<nFeatures;i++)
        {
            if(features[i-1].stop<bedPEEntry.start2 && bedPEEntry.start2<features[i].start)
            {
                strandCountPair.count+=features[i].startCount; return;
            }
        }
    });

    string name2;
    string strand2;

    if(countByGene2.empty()==false)
    {
        auto entry=countByGene2.begin();
        auto entryEnd=countByGene2.end();

        for(;entry!=entryEnd;entry++)
        {
            if(gotExonic2==true && entry->second.isExonic==false) continue;
            name2=entry->first+string(intronExon[entry->second.isExonic]);
            strand2=entry->second.strand;
            break;
        }

        for(entry++;entry!=entryEnd;entry++)
        {
            if(gotExonic2==true && entry->second.isExonic==false) continue;
            name2+=string(",")+entry->first+string(intronExon[entry->second.isExonic]);
            strand2+=string(",")+entry->second.strand;
        }
    }
    else
    {
        auto & boundariesByTarget=geneBoundaries[bedPEEntry.tid2];
        auto downstreamGene=boundariesByTarget.upper_bound(bedPEEntry.start2);

        if(downstreamGene!=boundariesByTarget.end())
        {
            int32_t disToDownstream=downstreamGene->first-bedPEEntry.start2;

            if(downstreamGene!=boundariesByTarget.begin())
            {
               auto upstreamGene=downstreamGene;--upstreamGene;

               int32_t disToUpstream=bedPEEntry.start2-upstreamGene->first;

                if(disToUpstream<disToDownstream)
                {
                    name2=upstreamGene->second+string("(")+to_string(disToUpstream)+string(")");
                }
                else
                {
                    name2=string("(")+to_string(disToDownstream)+string(")")+downstreamGene->second;
                }
            }
            else
            {
                name2=string("(")+to_string(disToDownstream)+string(")")+downstreamGene->second;
            }
        }
        else
        {
            if(downstreamGene!=boundariesByTarget.begin())
            {
                auto upstreamGene=downstreamGene;--upstreamGene;
                name2=upstreamGene->second+string("(")+to_string(bedPEEntry.start2-upstreamGene->first)+string(")");
            }
        }

        strand2=string(".");
    }

    //----------------------------------------------------------------
    //Compute score
    //----------------------------------------------------------------

    double score=0.0;

    if(bedPEEntry.splitCount==0) goto WRITE_ENTRY;

    if(gotExonic1==true && gotExonic2==false)
    {
        int32_t posCount1=0;
        int32_t negCount1=0;
        for(const auto & entry : countByGene1)
        {
            if(entry.second.strand=='+') posCount1+=entry.second.count;
            if(entry.second.strand=='-') negCount1+=entry.second.count;
        }

        if(posCount1){ score=double(bedPEEntry.splitCount)/double(bedPEEntry.splitCount+posCount1); goto WRITE_ENTRY;}
        if(negCount1){ score=double(bedPEEntry.splitCount)/double(bedPEEntry.splitCount+negCount1); goto WRITE_ENTRY;}

        int32_t posCount2=0;
        int32_t negCount2=0;
        for(const auto & entry : countByGene2)
        {
            if(entry.second.strand=='+') posCount2+=entry.second.count;
            if(entry.second.strand=='-') negCount2+=entry.second.count;
        }

        if(negCount2){ score=double(bedPEEntry.splitCount)/double(bedPEEntry.splitCount+negCount2); goto WRITE_ENTRY;}
        score=double(bedPEEntry.splitCount)/double(bedPEEntry.splitCount+posCount2); goto WRITE_ENTRY;
    }

    if(gotExonic1==false && gotExonic2==true)
    {
        int32_t posCount2=0;
        int32_t negCount2=0;
        for(const auto & entry : countByGene2)
        {
            if(entry.second.strand=='+') posCount2+=entry.second.count;
            if(entry.second.strand=='-') negCount2+=entry.second.count;
        }

        if(negCount2){ score=double(bedPEEntry.splitCount)/double(bedPEEntry.splitCount+negCount2); goto WRITE_ENTRY;}
        if(posCount2){ score=double(bedPEEntry.splitCount)/double(bedPEEntry.splitCount+posCount2); goto WRITE_ENTRY;}

        int32_t posCount1=0;
        int32_t negCount1=0;
        for(const auto & entry : countByGene1)
        {
            if(entry.second.strand=='+') posCount1+=entry.second.count;
            if(entry.second.strand=='-') negCount1+=entry.second.count;
        }

        if(posCount1){ score=double(bedPEEntry.splitCount)/double(bedPEEntry.splitCount+posCount1); goto WRITE_ENTRY;}
        score=double(bedPEEntry.splitCount)/double(bedPEEntry.splitCount+negCount1); goto WRITE_ENTRY;
    }

    if(gotExonic1==gotExonic2)
    {
        int32_t posCount1=0;
        int32_t negCount1=0;
        for(const auto & entry : countByGene1)
        {
            if(entry.second.strand=='+') posCount1+=entry.second.count;
            if(entry.second.strand=='-') negCount1+=entry.second.count;
        }

        int32_t posCount2=0;
        int32_t negCount2=0;
        for(const auto & entry : countByGene2)
        {
            if(entry.second.strand=='+') posCount2+=entry.second.count;
            if(entry.second.strand=='-') negCount2+=entry.second.count;
        }

        if(posCount1){ score=double(bedPEEntry.splitCount)/double(bedPEEntry.splitCount+posCount1); goto WRITE_ENTRY;}
        if(negCount2){ score=double(bedPEEntry.splitCount)/double(bedPEEntry.splitCount+negCount2); goto WRITE_ENTRY;}
        if(negCount1){ score=double(bedPEEntry.splitCount)/double(bedPEEntry.splitCount+negCount1); goto WRITE_ENTRY;}
        score=double(bedPEEntry.splitCount)/double(bedPEEntry.splitCount+posCount2);
    }

    //----------------------------------------------------------------
    //Write bedPE entry
    //----------------------------------------------------------------

WRITE_ENTRY:

    int32_t distJunction=-1; if(bedPEEntry.tid1==bedPEEntry.tid2) distJunction=bedPEEntry.start2-bedPEEntry.stop1;

    bedPEFile << bamHeader->target_name[bedPEEntry.tid1] << '\t'
              << bedPEEntry.start1 << '\t'
              << bedPEEntry.stop1  << '\t'
              << bamHeader->target_name[bedPEEntry.tid2] << '\t'
              << bedPEEntry.start2 << '\t'
              << bedPEEntry.stop2 << '\t'
              << name1+" - "+name2 << '\t'
              << score << '\t'
              << strand1 << '\t'
              << strand2 << '\t'
              << distJunction << '\t'
              << bedPEEntry.splitCount << '\t'
              << bedPEEntry.spanningCount << '\n';
}
//----------------------------------------------------------------
bool AnomalousJunctionCount::WriteBedPEFile(const bam_hdr_t * bamHeader)
{
    //----------------------------------------------------------------
    //Create junctions file
    //----------------------------------------------------------------

    ofstream bedPEFile(outputPrefix+string("_junctions.bedpe"),ofstream::out);

    if(bedPEFile.good()==false)
    {
        cerr << "Error: Could not create junctions bedpe file!" << endl;
        return false;
    }

    bedPEFile << "chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tname\tscore\tstrand1\tstrand2\tjunctionDist\tsplit\tspanning\n";

    //----------------------------------------------------------------
    //Iterate over split junctions
    //----------------------------------------------------------------

    BedPEEntry bedPEEntry;

    int32_t nTargets=bamHeader->n_targets;

    for(int32_t tid=0;tid<nTargets;tid++)
    {
        bedPEEntry.tid1=bedPEEntry.tid2=tid;

        for(const auto & juncPos1 : splitJunctions[tid])
        {
            for(const auto & juncPos2 : juncPos1.second)
            {
                if(juncPos2.second.splitCount+juncPos2.second.spanningCount<minJunctionCount) continue;

                if(juncPos2.second.splitCount==0)
                {
                    bool inSingleExon=false;

                    intervalTree[tid].visit_overlapping(juncPos1.first,juncPos2.first,[&](const Interval<int32_t,IntervalTreeEntry> & entry)
                    {
                        if(entry.start<=juncPos1.first && juncPos2.first<=entry.stop)
                        {
                            inSingleExon=true;
                            return;
                        }
                    });

                    if(inSingleExon) continue; //Skip if in single exon and spanning only (Extra filter from mathijs)
                }

                bedPEEntry.start1=juncPos2.second.start1;
                bedPEEntry.stop1=juncPos1.first;
                bedPEEntry.start2=juncPos2.first;
                bedPEEntry.stop2=juncPos2.second.stop2;
                bedPEEntry.splitCount=juncPos2.second.splitCount;
                bedPEEntry.spanningCount=juncPos2.second.spanningCount;

                WriteBedPEEntry(bamHeader,bedPEFile,bedPEEntry);
            }
        }
    }

    //----------------------------------------------------------------
    //Write interchromosomal junctions
    //----------------------------------------------------------------

    for(const auto & juncPos1 : interchromosomalJunctions)
    {
        bedPEEntry.tid1=juncPos1.first;

        for(const auto & juncPos2 : juncPos1.second)
        {
            bedPEEntry.tid2=juncPos2.first;

            for(const auto & junction : juncPos2.second)
            {
                if(junction.count<minJunctionCount) continue;

                bedPEEntry.start1=junction.start1;
                bedPEEntry.stop1=junction.stop1;
                bedPEEntry.start2=junction.start2;
                bedPEEntry.stop2=junction.stop2;
                bedPEEntry.splitCount=0;
                bedPEEntry.spanningCount=junction.count;

                WriteBedPEEntry(bamHeader,bedPEFile,bedPEEntry);
            }
        }
    }

    if(bedPEFile.good()==false)
    {
        cerr << "Error: Could not write bedpe file!" << endl;
        return false;
    }

    //----------------------------------------------------------------
    //Done
    //----------------------------------------------------------------

    return true;
}
//----------------------------------------------------------------
bool AnomalousJunctionCount::SecondPass(void)
{
    //----------------------------------------------------------------
    //Open spanning bam file
    //----------------------------------------------------------------

    cerr << "Info: Open spanning bam file" << endl;

    string spanningBamFilename=outputPrefix+"_spanning.bam";
    htsFile * spanningBamFile=sam_open(spanningBamFilename.c_str(),"r");

    if(spanningBamFile==nullptr)
    {
        cerr << "Error: Could not open bam file: " << spanningBamFilename << endl;
        return false;
    }

    //----------------------------------------------------------------
    //Read header spanning bam file
    //----------------------------------------------------------------

    bam_hdr_t * bamHeader=sam_hdr_read(spanningBamFile);

    if(bamHeader==nullptr)
    {
        cerr << "Error: Could not read bam header: " << spanningBamFilename << endl;
        sam_close(spanningBamFile); return false;

    }

    bam1_t * read1=nullptr;
    bam1_t * read2=bam_init1();

    //----------------------------------------------------------------
    //Initialize spanning junctions
    //----------------------------------------------------------------

    cerr << "Info: Initialize spanning junctions" << endl;

    size_t nTargets=bamHeader->n_targets;
    vector<multimap<int32_t,SpanningJunction> > spanningJunctions(nTargets);

    //----------------------------------------------------------------
    //Process spanning bam file
    //----------------------------------------------------------------

    cerr << "Info: Process spanning bam file" << endl;

    bool ret=true;

    string queryName1;
    string queryName2;
    vector<Segment> segments;

    while(sam_read1(spanningBamFile,bamHeader,read2)>=0)
    {
        //----------------------------------------------------------------
        //Group reads. Reads are already grouped by name
        //----------------------------------------------------------------

        if(read1==nullptr){read1=bam_dup1(read2);continue;}

        queryName1=string(bam_get_qname(read1))+':'+to_string(read1->core.pos)+':'+to_string(read1->core.mpos);
        queryName2=string(bam_get_qname(read2))+':'+to_string(read2->core.mpos)+':'+to_string(read2->core.pos);

        if(queryName1!=queryName2)  //Reads should be ordered by name in spanning bam file!
        {
            cerr << "Error: Reads in spanning bam file should be ordered by name!" << endl;
            ret=false; goto CLEAN_UP;
        }

        int32_t tid=read1->core.tid;
        if(tid!=read2->core.tid)
        {
            cerr << "Error: Target IDs of read mates should agree in the spanning bam file!" << endl;
            ret=false; goto CLEAN_UP;
        }

        //----------------------------------------------------------------
        //Create segments
        //----------------------------------------------------------------

        size_t nSegments1,nSegments2;
        if(CreateSegments(read1,read2,nSegments1,nSegments2,segments)==false)
        {
            cerr << "Error: Reads in spanning bam file should not be inconsistent!" << endl;
            ret=false; goto CLEAN_UP;
        }

        size_t nSegments=segments.size();
        if(nSegments1+nSegments2!=nSegments)
        {
            cerr << "Error: Reads in spanning bam file should never overlap!" << endl;
            ret=false; goto CLEAN_UP;
        }

        int32_t initFragmentLen=0; for(size_t i=0;i<nSegments;i++) initFragmentLen+=1 + segments[i].stop - segments[i].start;  //Calculate initial fragment length (This fragment is at least this long)

        //----------------------------------------------------------------
        //Add spanning junction to split junctions
        //----------------------------------------------------------------

        const auto & segment1=segments[nSegments1-1];
        const auto & segment2=segments[nSegments1];
        SpanningJunction spanningJunction(segment1.start,segment1.stop,segment2.start,segment2.stop);

        if(AddSpanToSplitJunctions(initFragmentLen,spanningJunction,splitJunctions[tid])==true)
        {
            bam_destroy1(read1);read1=nullptr;
            continue;
        }

        //----------------------------------------------------------------
        //Merge spanning junctions
        //----------------------------------------------------------------

        MergeSpanningJunctions(spanningJunction,spanningJunctions[tid]);
        bam_destroy1(read1); read1=nullptr;
    }

    //----------------------------------------------------------------
    //Add spanning junctions to split junctions
    //----------------------------------------------------------------

    cerr << "Info: Add spanning junctions to split junctions" << endl;

    for(size_t tid=0;tid<nTargets;tid++)
    {
        for(const auto & itSpanningJunction : spanningJunctions[tid])
        {
            const auto spanningJunction=itSpanningJunction.second;
            splitJunctions[tid][spanningJunction.stop1][spanningJunction.start2].Init(spanningJunction);
        }
    }

    //----------------------------------------------------------------
    //Write bedPE file
    //----------------------------------------------------------------

    cerr << "Info: Write bedPE file" << endl;

    ret=WriteBedPEFile(bamHeader);

    //----------------------------------------------------------------
    //CLean up
    //----------------------------------------------------------------

CLEAN_UP:

    if(read1!=nullptr)bam_destroy1(read1);
    bam_destroy1(read2);
    sam_hdr_destroy(bamHeader);
    sam_close(spanningBamFile);

    remove((outputPrefix+"_spanning.bam").c_str()); //Remove spanning bam file

    //----------------------------------------------------------------
    //Done
    //----------------------------------------------------------------

    cerr << "Info: Done" << endl;

    return ret;
}
//----------------------------------------------------------------
AnomalousJunctionCount::AnomalousJunctionCount(void) : countDuplicates(false),countSupplementary(false),maxFragLen(500),minJunctionCount(6),minMappingQuality(40),transcriptExtension(20){}
//----------------------------------------------------------------
#define INPUT_BAM_FILE          'b'
#define GTF_FILE                'g'
#define OUTPUT_PREFIX           'o'
#define MAX_FRAG_LEN            'l'
#define MIN_JUNCTION_COUNT      'c'
#define MIN_MAPPING_QUALITY     'q'
#define TRANSCRIPT_EXTENSION    'e'
#define COUNT_DUPLICATES        'D'
#define COUNT_SUPPLEMENTARY     's'
#define HELP                    'h'
#define SHORT_OPTIONS           "b:g:o:l:c:q:e:Dsh"
//----------------------------------------------------------------
struct option longOptions[] =
{
    {"input-bam-file",required_argument,nullptr,INPUT_BAM_FILE},
    {"gtf-file",required_argument,nullptr,GTF_FILE},
    {"output-prefix",required_argument,nullptr,OUTPUT_PREFIX},
    {"max-frag-len",required_argument,nullptr,MAX_FRAG_LEN},
    {"min-junction-count",required_argument,nullptr,MIN_JUNCTION_COUNT},
    {"min-mapping-quality",required_argument,nullptr,MIN_MAPPING_QUALITY},
    {"transcript-extension",required_argument,nullptr,TRANSCRIPT_EXTENSION},
    {"count-duplicates",no_argument,nullptr,COUNT_DUPLICATES},
    {"count-supplementary",no_argument,nullptr,COUNT_SUPPLEMENTARY},
    {"help",no_argument,nullptr,HELP}
};
//----------------------------------------------------------------
int AnomalousJunctionCount::Run(int argc,char * argv[])
{
    //----------------------------------------------------------------
    //Get input arguments
    //----------------------------------------------------------------

    cerr << "Info: Get input arguments" << endl;

    bool showHelp=(argc==1);
    int option,optionIndex;

    while((option=getopt_long(argc,argv,SHORT_OPTIONS,longOptions,&optionIndex))>=0)
    {
        switch(option)
        {
            case INPUT_BAM_FILE: inputBamFilename=string(optarg); continue;
            case GTF_FILE: gtfFilename=string(optarg); continue;
            case OUTPUT_PREFIX: outputPrefix=string(optarg);continue;
            case MAX_FRAG_LEN: maxFragLen=atoi(optarg); continue;
            case MIN_JUNCTION_COUNT: minJunctionCount=atoi(optarg); continue;
            case MIN_MAPPING_QUALITY: minMappingQuality=uint8_t(min(max(atoi(optarg),0),255)); continue;
            case TRANSCRIPT_EXTENSION: transcriptExtension=atoi(optarg); continue;
            case COUNT_DUPLICATES: countDuplicates=true; continue;
            case COUNT_SUPPLEMENTARY: countSupplementary=true; continue;
            case HELP: showHelp=true; continue;
        }
    }

    //----------------------------------------------------------------
    //Show help
    //----------------------------------------------------------------

    if(showHelp)
    {
        cerr << "anomalous_junction_count [options]"                                                        << endl;
        cerr                                                                                                << endl;
        cerr << "-b --input-bam-file <text>         Single input bam file (required)"                       << endl;
        cerr << "-g --gtf-file <text>               Single GTF file (required)"                             << endl;
        cerr << "-o --output-prefix <text>          Output prefix (optional default bam without extension)" << endl;
        cerr << "-l --max-frag-len <int>            Maximum fragment length (optional default = 500)"       << endl;
        cerr << "-c --min-junction-count <int>      Minimum junction count (optional default = 10)"         << endl;
        cerr << "-q --min-mapping-quality <int>     Minimum mapping quality (optional default = 40)"        << endl;
        cerr << "-e --transcript-extension <int>    Transcript extension (optional default = 20)"           << endl;
        cerr << "-D --count-duplicates <void>       Count duplicates when specified"                        << endl;
        cerr << "-s --count-supplementary <void>    Count supplementary when specified"                     << endl;
        cerr << "-h --help <void>                   This help"                                              << endl;
        cerr                                                                                                << endl;

        return 0;
    }

    //----------------------------------------------------------------
    //Check input arguments
    //----------------------------------------------------------------

    cerr << "Info: Check input arguments" << endl;

    if(inputBamFilename.empty())
    {
        cerr << "Error: Please specify a bam file" << endl;
        return 1;
    }

    if(gtfFilename.empty())
    {
        cerr << "Error: Please specify a gtf file" << endl;
        return 1;
    }

    if(outputPrefix.empty()) outputPrefix=inputBamFilename.substr(0,inputBamFilename.find_last_of('.'));    //Output prefix is input bam file without the extension
    if(transcriptExtension<0) transcriptExtension=0;

    //----------------------------------------------------------------
    //Run first pass (Find known, split and interchromosomal junctions)
    //----------------------------------------------------------------

    if(FirstPass()==false) return 1;

    //----------------------------------------------------------------
    //Run second pass (Add spanning junctions and write bedPE file)
    //----------------------------------------------------------------

    if(SecondPass()==false) return 1;

    //----------------------------------------------------------------
    //Done
    //----------------------------------------------------------------

    return 0;
}
//----------------------------------------------------------------
/*
if(string(bam_get_qname(read2))=="A00383:749:HJVL2DRX5:1:2145:9353:34068")
{
    cout << "Howdy" << endl;
}
*/
