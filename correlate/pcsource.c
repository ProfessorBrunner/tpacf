#include "correlation.h"

/************************************************************************************

 Contains routines to allocate space for, read in, and free angular and spatial data

 ************************************************************************************/

/***** Utility routines for parallel and serial code *****/

void pcsourceGetMemory(pcsource **block1,pcsource **block2,char *fileList){

/* Allocates two blocks of memory large enough to hold the actual xyz data */

	char sourceName[BUFFER_SIZE],fileName[BUFFER_SIZE];
	int i,NumSamples,NumData,NumFiles,max;
	FILE *list,*source;
	
	max = 0;
	
	#ifdef USE_MPI
	if(MyRank == 0){
	#endif
		list = fopen(fileList,"r");
		if(list == NULL){
			fprintf(stderr,"Failed to open file %s\n",fileList);
			return;
		}
		while(!feof(list)){
			fscanf(list,"%s %d",sourceName,&NumFiles);
			if(!feof(list)){
				for(i=0;i<NumFiles;i++){
					sprintf(fileName,"%s%d.bin",sourceName,i);
					source = fopen(fileName,"r");
					if(source == NULL){
						fprintf(stderr,"Failed to open file %s\n",fileName);
						return;
					}
					fread(&NumSamples,sizeof(int),1,source);
					fread(&NumData,sizeof(int),1,source);
					fclose(source);
					if(NumData > max){
						max = NumData;
					}
				}
			}
		}
		fclose(list);
		
	#ifdef USE_MPI
		MPI_Bcast(&max,1,MPI_INT,0,MPI_COMM_WORLD);
	}else{
		MPI_Bcast(&max,1,MPI_INT,0,MPI_COMM_WORLD);
	}
	#endif
	
	*block1 = malloc(max*sizeof(pcsource));
	*block2 = malloc(max*sizeof(pcsource));
	
	if(*block1 == NULL || *block2 == NULL){
		fprintf(stderr,"Failed to allocate memory for data/rands\n");
		exit(1);
	}
	
	return;
}


/******* Serial-only functions *******/

#ifndef USE_MPI
pcsource *pcsourceRead(char dataName[],pcsource *data,int **Samples){

/* Reads in points from dataName.bin and returns pointer to them.	*
 * Assumes data points to previously allocated space			*/

	char dataFile[BUFFER_SIZE];
	FILE *in;

	TIMESTART(startTime);
	
	/* Append the data sets name */
	sprintf(dataFile,"%s.bin",dataName);

	/* Open the file */
	if((in = fopen(dataFile,"r")) == NULL){
		fprintf(stderr,"Failed to open file %s\n",dataFile);
		return NULL;
	}
	
	/* Read the number of jackknife samples */
	if(fread(&NumSamples,sizeof(int),1,in) != 1){
		fprintf(stderr,"Failed to read header from %s\n",dataFile);
		exit(1);
	}

	/* Allocate space to store sample counts, if necessary */
	if(*Samples == NULL){
		if((*Samples=malloc((NumSamples+1)*sizeof(int))) == NULL){
			fprintf(stderr,"Failed to allocate sample space\n");
			exit(1);
		}
	}
	
	/* Read in sample counts */
	if(fread(*Samples,sizeof(int),(NumSamples+1),in) != (NumSamples+1)){
		fprintf(stderr,"Failed to read header from %s\n",dataFile);
		exit(1);
	}
	
	/* Check if data points to a valid memory address (sort of) */
	if(data == NULL){
		fprintf(stderr,"pcsource memory corrupted...aborting\n");
		exit(1);
	}
	
	/* Finally read in the actual data */
	if(fread(data,sizeof(pcsource),**Samples,in) != **Samples){
		fprintf(stderr,"Failed to read data from %s\n",dataFile);
		exit(1);
	}
	
	fclose(in);

	TIMESTOP(wallDataRead,startTime);
	
	return data;
}
#endif


/****** MPI-only functions ******/

#ifdef USE_MPI
static void pcsourceBuildMPI(MPI_Datatype *MPI_pcsource){

/* Builds MPI derived type corresponding to C struct pcsource */

	int BlockLengths[3];
	MPI_Datatype typelist[3];
	MPI_Aint StartAddress,Address;
	MPI_Aint Displacements[3];
	pcsource protopcsource;

	BlockLengths[0] = BlockLengths[1] = BlockLengths[2] = 1;
	typelist[0] = typelist[1] = typelist[2] = MPI_DOUBLE;

	Displacements[0] = 0;
	MPI_Address(&(protopcsource.x),&StartAddress);
	MPI_Address(&(protopcsource.y),&Address);
	Displacements[1] = Address - StartAddress;
	MPI_Address(&(protopcsource.z),&Address);
	Displacements[2] = Address - StartAddress;

	MPI_Type_struct(3,BlockLengths,Displacements,typelist,MPI_pcsource);
	MPI_Type_commit(MPI_pcsource);

	return;// MPI_pcsource;
}

pcsource *pcsourceRead(char dataName[],pcsource *data,int **Samples){

	static int firstc = 1;
	static MPI_Datatype MPI_pcsource;
	char dataFile[BUFFER_SIZE];
	int ierr,NumRead;
	MPI_File in;
	MPI_Status status;
	
	TIMESTART(startTime);
	
	if(firstc == 1){
		pcsourceBuildMPI(&MPI_pcsource);
		firstc = 0;
	}

	DEBUGPRINT("Entering pcsourceRead\n");

	sprintf(dataFile,"%s.bin",dataName);
	/* Open the file */
	if((ierr=MPI_File_open(MPI_COMM_WORLD,&dataFile[0],MPI_MODE_RDONLY,MPI_INFO_NULL,&in)) != MPI_SUCCESS){
		fprintf(stderr,"Failed to open data file %s\n",dataFile);
		exit(1);
	}
	
	/* Get the number of jackknife samples */
	if((ierr=MPI_File_read_all(in,&NumSamples,1,MPI_INT,&status)) != MPI_SUCCESS){
		fprintf(stderr,"Failed to read header with error %d\n",ierr);
		exit(1);
	}
	/* Allocate space to store sample counts, if necessary */
	if(*Samples == NULL){
		if((*Samples=malloc((NumSamples+1)*sizeof(int))) == NULL){
			fprintf(stderr,"Failed to allocate sample space\n");
			exit(1);
		}
	}
	/* Now read in sample counts */
	if((ierr=MPI_File_read_all(in,*Samples,(NumSamples+1),MPI_INT,&status)) != MPI_SUCCESS){
		fprintf(stderr,"Failed to read header with error %d\n",ierr);
		exit(1);
	}

	/* Read the data */
	if((ierr=MPI_File_read_all(in,&data[0],**Samples,MPI_pcsource,&status)) != MPI_SUCCESS){
		fprintf(stderr,"Failed to read data with error %d\n",ierr);
		exit(1);
	}

	/* Make sure we got it all */
	MPI_Get_count(&status,MPI_pcsource,&NumRead);
	if(NumRead != **Samples){
		fprintf(stderr,"MPI_File_read only got %d of %d data points\n",NumRead,**Samples);	
		exit(1);
	}

	MPI_File_close(&in);
	
	DEBUGPRINT("Exiting pcsourceRead\n");
	
	TIMESTOP(wallDataRead,startTime);
	
	return data;
}
#endif
