//----------------------------------------------------------------
// Name        : main.cpp
// Author      : Remco Hoogenboezem
// Version     :
// Copyright   :
// Description :
//----------------------------------------------------------------
#include "anomalous_junction_count.h"
//----------------------------------------------------------------
int main(void)
{
    char command[]="/scripts/rnapipe/bin/anomalous_junction_count6 -q 1 -b 12732.bam -g /data/valk_group/aszabo/RNA/reanalysis/anomalous_junctions/tools/combined.gtf.gz -o 12732_noMAPQ -D -c 1";
    int argc=12; char * argv[12]; char * pCommand=command; for(int i=0;i<argc;i++) argv[i]=strsep(&pCommand," ");

    AnomalousJunctionCount anomalousJunctionCount;
    return anomalousJunctionCount.Run(argc,argv);
}
//----------------------------------------------------------------


