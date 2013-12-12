#include "correlation.h"

bin *binBuildSpatial(double minDist,double maxDist){

/* Allocate bins and initialize their parameters.  Bins are 			*
 * logarithmically spaced in distance, but binning is done in distance^2	*
 *										*
 * RETURN:  pointer to the array of bins					*/	

	int i,j;
	bin *bins;
	
	DEBUGPRINT("Entering binBuildSpatial\n");

	bins = malloc((NumBins+2)*sizeof(bin));
	if(bins == NULL){
		printf("Failed to allocate memory in binBuild\n");
		return NULL;
	}

	for(i=0;i<=NumBins;i++){
		bins[i].limit  = maxDist*maxDist*pow(minDist/maxDist,(2.*(double)i)/((double)NumBins));
		bins[i].center  = maxDist*pow(minDist/maxDist,((double)i-0.5)/((double)NumBins));
		bins[i].Cnt = malloc((NumSamples+1)*sizeof(unsigned long long int));
		if(bins[i].Cnt == NULL){
			fprintf(stderr,"Failed to allocate sample space in binBuild\n");
			return NULL;
		}
		for(j=0;j<=NumSamples;j++){
			bins[i].Cnt[j] = 0;
		}
	}
	bins[NumBins+1].limit = 0.;
	bins[NumBins+1].Cnt = malloc((NumSamples+1)*sizeof(unsigned long long int));
	if(bins[NumBins+1].Cnt == NULL){
		fprintf(stderr,"Failed to allocate sample space in binBuild\n");
		return NULL;
	}
	for(i=0;i<=NumSamples;i++){
		bins[NumBins+1].Cnt[i] = 0;
	}
	
	DEBUGPRINT("Exiting binBuildSpatial\n");

	return bins;
}

bin *binBuildAngular(double minAngle,double maxAngle){

/* Allocate bins and initialize their parameters.  Bins are 			*
 * logarithmically spaced in theta, but binning is done in cos(theta)		*
 *										*
 * RETURN:  pointer to the array of bins					*/

	int i,j;
	bin *bins;

	DEBUGPRINT("Entering binBuildAngular\n");
	
	bins = malloc((NumBins+2)*sizeof(bin));
	if(bins == NULL){
		printf("Failed to allocate memory in buildbins\n");
		return NULL;
	}

	for(i=0;i<=NumBins;i++){
		bins[i].limit  = cos(maxAngle*pow(minAngle/maxAngle,((double)i)/((double)NumBins)));
		/* Get the center of the bin in degrees for printing */
		bins[i].center  = 180./M_PI*maxAngle*pow(minAngle/maxAngle,((double)i-0.5)/((double)NumBins));
		bins[i].Cnt = malloc((NumSamples+1)*sizeof(unsigned long long int));
		if(bins[i].Cnt == NULL){
			fprintf(stderr,"Failed to allocate sample space in binBuild\n");
			return NULL;
		}
		for(j=0;j<=NumSamples;j++){
			bins[i].Cnt[j] = 0;
		}
	}
	bins[NumBins+1].limit = 2.;
	bins[NumBins+1].Cnt = malloc((NumSamples+1)*sizeof(unsigned long long int));
	if(bins[NumBins+1].Cnt == NULL){
		fprintf(stderr,"Failed to allocate sample space in binBuild\n");
		return NULL;
	}
	for(i=0;i<=NumSamples;i++){
		bins[NumBins+1].Cnt[i] = 0;
	}

	DEBUGPRINT("Exiting binBuildAngular\n");
	
	return bins;
}


void binPrint(char filename[],bin bins[],int Samples1[],int Samples2[],int BinType){

/*  Prints counts in bins to filename				*/

	int i,j;
	unsigned long long int MultFactor;
	FILE *out;

	out = fopen(filename,"w");
	if(out == NULL) printf("Failed to open %s in binPrint\n",filename);
	if(BinType == 0){
		MultFactor = 2;
		fprintf(out,"%d %d %d",NumBins,NumSamples,Samples1[0]);
		for(i=1;i<=NumSamples;i++){
			fprintf(out," %d",Samples1[i]);
			fflush(out);
		}
		fprintf(out,"\n");
	}else{
		MultFactor = 1;
		fprintf(out,"%d %d %d",NumBins,NumSamples,Samples1[0]);
		for(i=1;i<=NumSamples;i++){
			fprintf(out," %d",Samples1[i]);
		}
		fprintf(out,"\n%d",Samples2[0]);
		for(i=1;i<=NumSamples;i++){
			fprintf(out," %d",Samples2[i]);
		}
		fprintf(out,"\n");
	}

	for(i=1;i<=NumBins;i++){
		fprintf(out,"%16.12f %19llu",bins[i].center,MultFactor*bins[i].Cnt[0]);
		for(j=1;j<=NumSamples;j++){
			fprintf(out," %19llu",MultFactor*bins[i].Cnt[j]);
		}
		fprintf(out,"\n");
	}

	fclose(out);

	return;

}


void binClear(bin bins[]){

/* Reinitializes bin counts to zero, including for jackknife samples	*/

	int i,j;
	
	for(i=0;i<=NumBins+1;i++){
		for(j=0;j<=NumSamples;j++){
			bins[i].Cnt[j] = 0;
		}
	}
	
	return;
}


void binFree(bin *bins){

/*  Free space allocated for bins					*/

	int i;

	for(i=0;i<=NumBins+1;i++){
		free(bins[i].Cnt);
	}
	free(bins);

	return;
}

#ifdef USE_OMP
unsigned long long int *GlobalCnts;

void binReduceOmp(bin bins[]){

	int i,j,NumCnts,index,NSp1;

	NumCnts = (NumBins+1)*(NumSamples+1);
	NSp1 = NumSamples+1;

	#pragma omp master
	{
		GlobalCnts = malloc(NumCnts*sizeof(unsigned long long int));
		if(GlobalCnts == NULL){
			fprintf(stderr,"Failed to allocate GlobalCnts\n");
		}
		for(i=0;i<NumCnts;i++){
			GlobalCnts[i] = 0;
		}
	}
	#pragma omp barrier

	for(i=0;i<=NumBins;i++){
		for(j=0;j<=NumSamples;j++){
			index = i*NSp1 + j;
			#pragma omp atomic
				GlobalCnts[index] += bins[i].Cnt[j];
		}
	}

	#pragma omp barrier
	#pragma omp master
	{
	for(i=0;i<=NumBins;i++){
		for(j=0;j<=NumSamples;j++){
			index = i*NSp1 + j;
			bins[i].Cnt[j] = GlobalCnts[index];
		}
	}
	free(GlobalCnts);
	}

	return;
}
#endif
				



#ifdef USE_MPI

static void long_long_sum(unsigned long long int *invec,unsigned long long int *inoutvec,int *len,MPI_Datatype *type){

	int i;

	type = NULL;	/* Avoid compiler complaints about unused parameters */
	i = (int)type;	/* continued */
	
	for(i=0;i<(*len);i++){
		inoutvec[i] += invec[i];
	}

	return;
}


void binReduce(bin bins[]){

/*  Coordinates the combining of bin counts for distributed processes	*
 *									*
 *  In:	 Process's local array of bins					*
 *	 Number of bins							*
 *	 Number of samples						*
 *									*	
 *  Out: Combined bin counts in process 0's array of bins		*/

	int i,j,Nsp1,Nbp1,index;
	unsigned long long int *Cnts,*Bcnts;
	MPI_Op operator;

	TIMESTART(startTime);
	
	DEBUGPRINT("Entering binReduce\n");

	Nsp1 = NumSamples + 1;
	Nbp1 = NumBins + 1;
	
	Cnts = malloc((Nbp1*Nsp1)*sizeof(unsigned long long int));
	Bcnts = malloc((Nbp1*Nsp1)*sizeof(unsigned long long int)); 
	if(Cnts == NULL || Bcnts == NULL){
		fprintf(stderr,"Failed to allocate memory in binReduce\n");
		return;
	}
	
	for(i=0;i<Nbp1;i++){
		for(j=0;j<=NumSamples;j++){
			index = i*Nsp1 + j;
			Cnts[index] = 0;
			Bcnts[index] = bins[i].Cnt[j];
		}
	}

	MPI_Op_create((MPI_User_function *)long_long_sum,1,&operator);
	MPI_Reduce(&(Bcnts[0]),&(Cnts[0]),Nbp1*Nsp1,MPI_UNSIGNED_LONG_LONG,operator,0,MPI_COMM_WORLD);
	for(i=0;i<Nbp1;i++){
		for(j=0;j<=NumSamples;j++){
			bins[i].Cnt[j] = Cnts[i*Nsp1 + j];
		}
	}
	
	free(Cnts);
	free(Bcnts);

	DEBUGPRINT("Exiting binReduce\n");

	TIMESTOP(wallReduce,startTime);
	
	return;
}

#endif

	
