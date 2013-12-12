#include "estimate.h"


double ls(bin ddbins[],bin **drbins,bin **rrbins,int CurBin,int CurSamp,int NumFiles,int NumData[],int **NumRand){

	int i;
	double dd,dr,rr;
	
	//printf("%d %d %d %d\n",CurBin,CurSamp,ddbins[CurBin].Cnt[CurSamp],NumData[CurSamp]*NumData[CurSamp]);
	
	dr = 0.;
	rr = 0.;
	for(i=1;i<NumFiles;i++){
		dr += ((double)drbins[i][CurBin].Cnt[CurSamp])/((double)(((unsigned long long int)NumData[CurSamp])*((unsigned long long int)NumRand[i][CurSamp])));
		rr += ((double)rrbins[i][CurBin].Cnt[CurSamp])/((double)(((unsigned long long int)NumRand[i][CurSamp])*((unsigned long long int)NumRand[i][CurSamp])));
		//printf("%llu\t%d\t%d\t%g\n",drbins[i][CurBin].Cnt[CurSamp],NumData[CurSamp],NumRand[i][CurSamp],dr);
	}
	dr = dr/(NumFiles - 1.);
	rr = rr/(NumFiles - 1.);
	dd = ((double)ddbins[CurBin].Cnt[CurSamp])/((double)(((unsigned long long int)NumData[CurSamp])*((unsigned long long int)NumData[CurSamp])));
	
	//printf("Return from ls %f %f %f\n",dd,dr,rr);
	
	return (dd - 2.*dr + rr)/rr;
}

