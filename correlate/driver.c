#include "correlation.h"
#include "globals.h"

int main(int argc,char *argv[]){

/*
Driver routine to coordinate input/output and parallel correlation calculation

Inputs:		Parameter File

Outputs:	Unnormalized bin counts including jackknife resampling
*/

	/**************** Variable Declarations *********************/
	
	char fileList[BUFFER_SIZE],dataName[BUFFER_SIZE],randName[BUFFER_SIZE],binOutFile[BUFFER_SIZE],curFile[BUFFER_SIZE];
	int i,j,NumData,NumRand,WorkLevel,NumDataFiles,NumRandFiles,Terminate;
	int *Dsamples,*Rsamples,getDD,getDR,getRR;
	double minDist,maxDist;
	FILE *list;
	pcsource *data,*rand,*block1,*block2;
	void *rootDataTree,*rootRandTree,*tblock1,*tblock2;
	bin *ddbins,*drbins,*rrbins;
	
	/* parallel specific variables */
	#if defined (USE_MPI) || defined (USE_OMP)
		int NumWorkNodes;
		void **dworknodes,**rworknodes;
	#endif
	
	/* OpenMP specific variables */
	#ifdef USE_OMP
		int MyId;
	#endif
	
	/* Timing specific variables */
	#ifdef TIMING
		#ifdef USE_MPI
			double startTimeTotal;
		#else
			clock_t startTimeTotal;
		#endif
	#endif
	
	/******************** End Variable Declarations ******************/

	/* Initialize MPI */
	#ifdef USE_MPI
		#ifdef USE_OMP
			int ThreadSupport;
			MPI_Init_thread(&argc,&argv,MPI_THREAD_MULTIPLE,&ThreadSupport);
		#else
			MPI_Init(&argc,&argv);
		#endif
		MPI_Comm_rank(MPI_COMM_WORLD,&MyRank);
		MPI_Comm_size(MPI_COMM_WORLD,&NumProcs);
	#endif

	/* Check to see if a parameter file was specified */
	if(argc != 2){
		#ifdef USE_MPI
			if(MyRank == 0){
				fprintf(stderr,"The parameter file must be included as the only command line argument...ABORTING\n");
			}
			MPI_Finalize();
		#else
			fprintf(stderr,"The parameter file must be included as the only command line argument...ABORTING\n");
		#endif
		exit(1);
	}

	/********************************************************************************
	 First get all the parameters from argv[1] then allocate and initialize variables
	 ********************************************************************************/
	TIMESTART(startTimeTotal);

	/* Get parameters from command line supplied file */
	#ifdef USE_MPI
		getParams(argv[1],fileList,&getDD,&getDR,&getRR,&WorkLevel,&minDist,&maxDist);
	#else
		getParams(argv[1],fileList,&getDD,&getDR,&getRR,&WorkLevel,&minDist,&maxDist);
	#endif
	/* Initialize timing variables */
	#ifdef TIMING
		wallTimeAC = wallTimeCC = wallTimeTotal = 0.;
		wallDataRead = wallTreeRead = wallSetup = 0.;
		#ifdef USE_MPI
			wallReduce = wallComm = 0;
		#endif
		#ifdef USE_OMP
			wallReduceOmp = 0;
		#endif
	#endif


	#if defined (USE_MPI) || defined (USE_OMP)		/* Is there any parallelism? */
	/* Get the number of work nodes at WorkLevel, and allocate space for pointers to them */
	NumWorkNodes = 1 << WorkLevel;
	if(AngOrSpa == 0){
		dworknodes = malloc(NumWorkNodes*sizeof(angTreeNode *));
		rworknodes = malloc(NumWorkNodes*sizeof(angTreeNode *));
	}else{
		dworknodes = malloc(NumWorkNodes*sizeof(spaTreeNode *));
		rworknodes = malloc(NumWorkNodes*sizeof(spaTreeNode *));
	}
	if(dworknodes == NULL || rworknodes == NULL){
		fprintf(stderr,"Failed to allocate space for work nodes\n");
		exit(1);
	}
	#endif

	/* Allocate memory for data points up front */
	pcsourceGetMemory(&block1,&block2,fileList);
	/* Now allocate memory for trees */
	treeGetMemory(&tblock1,&tblock2,fileList);

	#ifdef USE_MPI
	if(MyRank == 0){
	#endif
		list = fopen(fileList,"r");
		fscanf(list,"%s %d",dataName,&NumDataFiles);
		sprintf(curFile,"%s%d",dataName,0);
	#ifdef USE_MPI
	}
	MPI_Bcast(&dataName[0],BUFFER_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&NumDataFiles,1,MPI_INT,0,MPI_COMM_WORLD);
	#endif

	/* Initialize sample counts to NULL */
	Dsamples = Rsamples = NULL;
	
	/* Begin OpenMP parallel region */
	#ifdef TIMING
	#pragma omp parallel private(ddbins,drbins,rrbins,MyId,i,j,startTime)
	#else
	#pragma omp parallel private(ddbins,drbins,rrbins,MyId,i,j)
	#endif
	{
	#ifdef USE_OMP
	MyId = omp_get_thread_num();
	#pragma omp master
	{ NumThreads = omp_get_num_threads(); }
	#endif

	#pragma omp master
	{ TIMESTOP(wallSetup,startTimeTotal); }
	
	/* Now start reading in data sets and getting requested counts */
	if(getDD == 1){
		#pragma omp master
		{ NumData = 0; }
		for(i=0;i<NumDataFiles;i++){
			#pragma omp master
			{
				sprintf(curFile,"%s%d",dataName,i);
				data = pcsourceRead(curFile,block1,&Dsamples);
				NumData += Dsamples[0];				// Add to total count for this data set
				rootDataTree = treeRead(curFile,tblock1);
				WORKNODES(rootDataTree,dworknodes);		// Fill list of worknodes, if parallel
			}
			#pragma omp barrier
			if(i==0) BINBUILD(ddbins);				// Build bins if necessary
			AC(data,dworknodes,rootDataTree,ddbins);		// Now count the pairs
			for(j=i+1;j<NumDataFiles;j++){				// Loop over this data sets files
				#pragma omp master
				{
					sprintf(curFile,"%s%d",dataName,j);
					rand = pcsourceRead(curFile,block2,&Rsamples);
					rootRandTree = treeRead(curFile,tblock2);
				}
				#pragma omp barrier
				CC(data,rand,dworknodes,rootDataTree,rootRandTree,ddbins);
			}
		}
		
		#ifdef USE_OMP
		binReduceOmp(ddbins);	// Sum up counts from shared memory processes
		#endif
		
		#ifdef USE_MPI
			#pragma omp master
			{ binReduce(ddbins); }	// Sum up counts from distributed processes
		#endif
		#pragma omp master
		{
			sprintf(binOutFile,"%s_ddbins",dataName);
			Dsamples[0] = NumData;
			#ifdef USE_MPI
			if(MyRank == 0)
			#endif
			binPrint(binOutFile,ddbins,Dsamples,Dsamples,0);	// Print out the results
		}
	}

	Terminate = -1;
	for( ; ; ){		/* Loop until all data sets have been processed */

		/* Get the next data name to processes */
		#pragma omp master
		{
			#ifdef USE_MPI
				if(MyRank == 0){
					fscanf(list,"%s %d",randName,&NumRandFiles);
					if(feof(list)) NumRandFiles = Terminate;
				}
				MPI_Bcast(&randName[0],BUFFER_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
				MPI_Bcast(&NumRandFiles,1,MPI_INT,0,MPI_COMM_WORLD);
			#else
				fscanf(list,"%s %d",randName,&NumRandFiles);
				if(feof(list)) NumRandFiles = Terminate;
			#endif
		}
		#pragma omp barrier
		
		/* Check if we're done */
		if(NumRandFiles == Terminate) break;
		
		#pragma omp master
		{
			sprintf(curFile,"%s0",randName);
			NumRand = 0;
		}
		
		/* Loop over all files for this data set */
		for(i=0;i<NumRandFiles;i++){
			#pragma omp master
			{
				sprintf(curFile,"%s%d",randName,i);
				rand = pcsourceRead(curFile,block1,&Rsamples);
				NumRand += Rsamples[0];
				rootRandTree = treeRead(curFile,tblock1);
				WORKNODES(rootRandTree,rworknodes);
			}
			if(i == 0){
				#pragma omp barrier
				BINBUILD(drbins);
			}
			if(getDR == 1){
			/* Get the DR counts */
				#pragma omp master
				{ NumData = 0; }
				for(j=0;j<NumDataFiles;j++){
					#pragma omp master
					{
						sprintf(curFile,"%s%d",dataName,j);
						data = pcsourceRead(curFile,block2,&Dsamples);
						NumData += Dsamples[0];
						rootDataTree = treeRead(curFile,tblock2);
					}
					#pragma omp barrier
					/* Get the counts */
					CC(rand,data,rworknodes,rootRandTree,rootDataTree,drbins);
				}
			}

			if(getRR == 1){
			/* Get the RR counts */
				#pragma omp barrier
				if(i == 0){
					BINBUILD(rrbins);
				}
				AC(rand,rworknodes,rootRandTree,rrbins);
				for(j=i+1;j<NumRandFiles;j++){
					#pragma omp master
					{
						sprintf(curFile,"%s%d",randName,j);
						data = pcsourceRead(curFile,block2,&Dsamples);
						rootDataTree = treeRead(curFile,tblock2);
					}
					#pragma omp barrier
					/* Get the counts */
					CC(rand,data,rworknodes,rootRandTree,rootDataTree,rrbins);
				}
			}
		}
			if(getDR == 1){
			/* Combine, print, and clear DR bin counts */
			#ifdef USE_OMP
				binReduceOmp(drbins);
			#endif
			#ifdef USE_MPI
				#pragma omp master
				{ binReduce(drbins); }
			#endif
				#pragma omp master
			{
				sprintf(binOutFile,"%s_%s_drbins",dataName,randName);
				Dsamples[0] = NumData;
				Rsamples[0] = NumRand;
				#ifdef USE_MPI
				if(MyRank == 0)
				#endif
				binPrint(binOutFile,drbins,Dsamples,Rsamples,1);
			}
			binClear(drbins);
		}
		if(getRR == 1){
			/* Combine, print, and clear RR bin counts */
			#ifdef USE_OMP
				binReduceOmp(rrbins);
			#endif
			#ifdef USE_MPI
				#pragma omp master
				{ binReduce(rrbins); }
			#endif
				#pragma omp master
				{
				sprintf(binOutFile,"%s_rrbins",randName);
				Rsamples[0] = NumRand;
				#ifdef USE_MPI
				if(MyRank == 0)
				#endif
				binPrint(binOutFile,rrbins,Rsamples,Rsamples,0);
				}
			binClear(rrbins);
		}
	}

	DEBUGPRINT("Cleaning up\n");
	
	if(getDD == 1){
		binFree(ddbins);
	}
	if(getDR == 1){
		binFree(drbins);
	}
	if(getRR == 1){
		binFree(rrbins);
	}

	}	/* End OpenMP parallel block */

	#ifdef USE_MPI
	if(MyRank == 0)
	#endif
		fclose(list);


	/* Free everything else */
	free(block1);
	free(block2);
	free(tblock1);
	free(tblock2);
	free(Dsamples);
	free(Rsamples);

	DEBUGPRINT("Reporting Timings\n");

	TIMESTOP(wallTimeTotal,startTimeTotal);
	/* Print the TIMING results */
	reportTimings();
	
	#ifdef USE_MPI
	/* Shut down MPI */	
	MPI_Finalize();
	#endif

	DEBUGPRINT("FINISHED!\n");
	
	/* Done! */
	return 0;
}

#ifdef USE_MPI
void getParams(char parameterFile[],char fileList[],int *getDD,int *getDR,int *getRR,int *WorkLevel,double *minDist,double *maxDist){

/*
Reads in parameters from parameterFile and distributes them among MPI processes

Note that this function converts the input angles from degrees to radians
*/

	//extern int NumBins;
	FILE *paramFile;
	
	if(MyRank == 0){
		paramFile = fopen(parameterFile,"r");

		fscanf(paramFile,"%s%*[^\n]",fileList);
		fscanf(paramFile,"%d%*[^\n]",getDD);
		fscanf(paramFile,"%d%*[^\n]",getDR);
		fscanf(paramFile,"%d%*[^\n]",getRR);
		fscanf(paramFile,"%d%*[^\n]",WorkLevel);
		fscanf(paramFile,"%d%*[^\n]",&NumBins);
		fscanf(paramFile,"%lf%*[^\n]",minDist);
		fscanf(paramFile,"%lf%*[^\n]",maxDist);
		fscanf(paramFile,"%d%*[^\n]",&AngOrSpa);

		if(AngOrSpa == 0){
			*minDist *= M_PI/180.;
			*maxDist *= M_PI/180.;
		}

		fclose(paramFile);
	}

	MPI_Bcast(&fileList[0],BUFFER_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(getDD,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(getDR,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(getRR,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(WorkLevel,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&NumBins,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(minDist,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(maxDist,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&AngOrSpa,1,MPI_INT,0,MPI_COMM_WORLD);

	return;
}
#else	/* NOT USE_MPI */
void getParams(char parameterFile[],char fileList[],int *getDD,int *getDR,int *getRR,int *WorkLevel,double *minDist,double *maxDist){

/*
Reads in parameters from parameterFile

Note that this function converts the input angles from degrees to radians
*/

	//extern int NumBins;
	FILE *paramFile;
	
	paramFile = fopen(parameterFile,"r");
	
	fscanf(paramFile,"%s%*[^\n]",fileList);
	fscanf(paramFile,"%d%*[^\n]",getDD);
	fscanf(paramFile,"%d%*[^\n]",getDR);
	fscanf(paramFile,"%d%*[^\n]",getRR);
	fscanf(paramFile,"%d%*[^\n]",WorkLevel);
	fscanf(paramFile,"%d%*[^\n]",&NumBins);
	fscanf(paramFile,"%lf%*[^\n]",minDist);
	fscanf(paramFile,"%lf%*[^\n]",maxDist);
	fscanf(paramFile,"%d%*[^\n]",&AngOrSpa);

	if(AngOrSpa == 0){
		*minDist *= M_PI/180.;
		*maxDist *= M_PI/180.;
	}
	
	fclose(paramFile);
	
	return;
}
#endif /* USE_MPI */

void reportTimings(){

	/*
	#if defined (USE_MPI) && defined (TIMING)
	double minRed,avgRed,maxRed,minComm,avgComm,maxComm;
		
	MPI_Reduce(&wallReduce,&minRed,1,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);
	MPI_Reduce(&wallReduce,&avgRed,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Reduce(&wallReduce,&maxRed,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
	MPI_Reduce(&wallComm,&minComm,1,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);
	MPI_Reduce(&wallComm,&avgComm,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Reduce(&wallComm,&maxComm,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
	avgRed = avgRed/NumProcs;
	avgComm = avgComm/NumProcs;
	#endif
	*/

	#ifdef USE_OMP
	wallTimeAC /= NumThreads;
	wallTimeCC /= NumThreads;
	#endif

	#ifdef TIMING
	#ifdef USE_MPI
		if(MyRank == 0){
			fprintf(stdout,"Total Walltime (s):\t%f\n\n",wallTimeTotal);
			fprintf(stdout,"AC - \ttime (s):\t%f\tperc. wall:\t%f\n",wallTimeAC,wallTimeAC/wallTimeTotal);
			fprintf(stdout,"CC - \ttime (s):\t%f\tperc. wall:\t%f\n",wallTimeCC,wallTimeCC/wallTimeTotal);
			fprintf(stdout,"Data Read - time (s):\t%f\tperc. wall:\t%f\n",wallDataRead,wallDataRead/wallTimeTotal);
			fprintf(stdout,"Tree Read - time (s):\t%f\tperc. wall:\t%f\n",wallTreeRead,wallTreeRead/wallTimeTotal);
			fprintf(stdout,"Setup - time (s):\t%f\tperc. wall:\t%f\n",wallSetup,wallSetup/wallTimeTotal);
			fprintf(stdout,"MPI Reduce - time (s):\t%f\tperc. wall:\t%f\n",wallReduce,wallReduce/wallTimeTotal);
			//fprintf(stdout,"\tavg - time (s):\t%f\tperc. wall:\t%f\n",avgRed,avgRed/wallTimeTotal);
			//fprintf(stdout,"\tmin - time (s):\t%f\tperc. wall:\t%f\n",minRed,minRed/wallTimeTotal);
			//fprintf(stdout,"\tmax - time (s):\t%f\tperc. wall:\t%f\n",maxRed,maxRed/wallTimeTotal);
			fprintf(stdout,"Work Req. - time(s):\t%f\tperc. wall:\t%f\n",wallComm,wallComm/wallTimeTotal);
			//fprintf(stdout,"\tavg - time (s):\t%f\tperc. wall:\t%f\n",avgComm,avgComm/wallTimeTotal);
			//fprintf(stdout,"\tmin - time (s):\t%f\tperc. wall:\t%f\n",minComm,minComm/wallTimeTotal);
			//fprintf(stdout,"\tmax - time (s):\t%f\tperc. wall:\t%f\n",maxComm,maxComm/wallTimeTotal);
		}
	#else
		fprintf(stdout,"Total Walltime (s):\t%f\n\n",wallTimeTotal);
		fprintf(stdout,"\tAC - %f\t\t%f\n",wallTimeAC,wallTimeAC/wallTimeTotal);
		fprintf(stdout,"\tCC - %f\t\t%f\n",wallTimeCC,wallTimeCC/wallTimeTotal);
		fprintf(stdout,"\tData Read - %f\t%f\n",wallDataRead,wallDataRead/wallTimeTotal);
		fprintf(stdout,"\tTree Read - %f\t%f\n",wallTreeRead,wallTreeRead/wallTimeTotal);	
	#endif

	#ifdef USE_OMP
	#ifdef USE_MPI
		if(MyRank == 0)
	#endif
		fprintf(stdout,"Omp Reduce - time (s):\t%f\tperc. wall:\t%f\n",wallReduceOmp,wallReduceOmp/wallTimeTotal);
	#endif
	#endif

	return;
}
