#include "precompute.h"

int main(int argc, char *argv[]){

/* Driver code to convert ra/dec pairs to x,y,z.  Also builds				*
 * tree for each data set contained in a file.  Converted data is output to a file(s) with	*
 * the same name as the original data file, appended with ".bin".  Tree is likewise	*
 * output to the same name, appended with ".tree".  Both are in unformatted binary,	*
 * appropriate for reading in the main angular correlation code.			*/

	int NumData,NumSamples,LogSamples,NumFiles,AngOrSpa;
	char infile[BUFFER_SIZE],outfile[BUFFER_SIZE];
	double dLogSamples;
	splitlog *dsplit;
	FILE *list,*newlist;
	
	if(argc != 4 && argc != 3){
		fprintf(stderr,"Usage:  precompute filenames.list 0 [NumSamples]\t for angular calculation\n");
		fprintf(stderr,"Usage:  precompute filenames.list 1 [NumSamples]\t for spatial calculation\n");
		exit(1);
	}

	AngOrSpa = atoi(argv[2]);

	if(AngOrSpa == 0){
		fprintf(stdout,"Performing angular precomputing\n");
	}else if(AngOrSpa == 1){
		fprintf(stdout,"Performing spatial precomputing\n");
	}else{
		fprintf(stdout,"Invalid argument\nAngular = 0\nSpatial = 1\n");
		exit(1);
	}
	
	if(argc == 4){
		NumSamples = atoi(argv[3]);
		if(NumSamples != 0){
			dLogSamples = log((double)NumSamples)/M_LN2;
			LogSamples = (int)(dLogSamples+0.5);
			if(fabs(dLogSamples-LogSamples) > ERR_TOL){
				printf("Number of Samples must be zero or a power of two...EXITING\n");
				exit(1);
			}
		}
	}else{
		NumSamples = DEFAULT_NUM_SAMPLES;
	}

	dsplit = malloc((NumSamples+1)*sizeof(splitlog));
	if(dsplit == NULL){
		printf("Failed to allocate splitlog...EXITING\n");
		exit(1);
	}

	list = fopen(argv[1],"r");
	if(list == NULL){
		fprintf(stderr,"Failed to open file list\n");
		exit(1);
	}
	fscanf(list,"%s%*[^\n]",&infile[0]);
	printf("%s ",infile);
	if(AngOrSpa == 0){
		NumData = convertAng(infile,NumSamples);
		NumFiles = angTree(infile,dsplit,NumData,NumSamples,0);
	}else{
		NumData = convertSpa(infile,NumSamples);
		NumFiles = spaTree(infile,dsplit,NumData,NumSamples,0);
	}
	if(NumFiles == 0){
		fprintf(stderr,"Failed to print results\n");
		exit(1);
	}
	sprintf(outfile,"%s.temp",argv[1]);
	newlist = fopen(outfile,"w");
	fprintf(newlist,"%s %d\n",infile,NumFiles);
	while(!feof(list)){	/* Loop until all files are processed */
		fscanf(list,"%s%*[^\n]",&infile[0]);
		if(!feof(list)){
			printf("%s ",infile);
			if(AngOrSpa == 0){
				NumData = convertAng(infile,NumSamples);
				NumFiles = angTree(infile,dsplit,NumData,NumSamples,1);
			}else{
				NumData = convertSpa(infile,NumSamples);
				NumFiles = spaTree(infile,dsplit,NumData,NumSamples,1);
			}
			if(NumFiles == 0){
				fprintf(stderr,"Failed to print results\n");
				exit(1);
			}
			fprintf(newlist,"%s %d\n",infile,NumFiles);
		}
	}
	free(dsplit);
	fclose(list);
	remove(argv[1]);
	fclose(newlist);
	rename(outfile,argv[1]);
	return 0;
}




