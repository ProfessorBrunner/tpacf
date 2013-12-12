#include "precompute.h"

/* Integration parameter */
#define REFINE_MAX 20

static double oneOverE(double z){

	double om,ol,zp1,zp1cubed;

	om = 0.3;	/* omega matter */
	ol = 0.7;	/* omega lambda */
	zp1 = z + 1.;
	zp1cubed = zp1*zp1*zp1;

	return 1./sqrt(om*zp1cubed+ol);
}	
		
static double romberg(double z){

	int i,j,steps,denom;
	double h,r1[REFINE_MAX+1],r2[REFINE_MAX+1],sum;
	
	h = z;
	r1[0] = oneOverE(0.);
	r1[0] += oneOverE(z);
	r1[0] *= 0.5*h;
	i = 1;
	steps = 1;
	denom = 4.;
	for( ; ; ){
		sum = 0.;
		for(j=1;j<=steps;j++){
			sum += oneOverE(((double)j-0.5)*h);
		}
		steps = steps << 1;
		r2[0] = .5*(r1[0] + h*sum);
		denom = 4;
		for(j=0;j<i;j++){
			r2[j+1] = r2[j] + (r2[j]-r1[j])/((double)denom-1.);
			denom = denom << 2;
		}
		if(fabs(r2[i] - r1[i-1]) < ERR_TOL*r1[i-1]){
			return r2[i];
		}else if(i == REFINE_MAX){
			fprintf(stderr,"Maximum refinement reached in romberg...ABORTING\n");
			exit(1);
		}
		h = h/2.;
		i++;
		for(j=0;j<i;j++){
			r1[j] = r2[j];
		}
	}
}
	
static double dm(double z){

	return 3000.*romberg(z);

}

int convertSpa(char infile[],int NumSamples){

/* Reads in data from file infile, converts it to x,y,z and outputs it to infile.bin.		*
 * Also prints information about the number of jackknife samples and number of data points.	*/

	int Cnt;
	char outfile[BUFFER_SIZE];
	double raRad,decRad,deg2rad,dist,cdecRad;
	zsource sdata;
	pcsource pcdata;
	FILE *in,*out;
	#ifndef USE_DISK
	int i;
	#endif
	
	in = fopen(infile,"r");
	#ifdef USE_DISK
	sprintf(outfile,"%s.bin.temp",infile);
	#else
	sprintf(outfile,"%s0.bin",infile);
	#endif
	out = fopen(outfile,"w");
	Cnt=NumSamples;
	
	#ifndef USE_DISK
	fwrite(&Cnt,sizeof(int),1,out);
	for(i=0;i<=NumSamples;i++){
		fwrite(&Cnt,sizeof(int),1,out);
	}
	#endif
	deg2rad = M_PI/180.;
	Cnt = 0;
	while(!feof(in)){
		fscanf(in,"%lf %lf %lf%*[^\n]",&sdata.ra,&sdata.dec,&sdata.z);
		if(!feof(in)){
			Cnt++;
			dist = dm(sdata.z);
			raRad = deg2rad*sdata.ra;
			decRad = deg2rad*sdata.dec;
			cdecRad = cos(decRad);
			pcdata.x = dist*cos(raRad)*cdecRad;
			pcdata.y = dist*sin(raRad)*cdecRad;
			pcdata.z = dist*sin(decRad);
			fwrite(&pcdata,sizeof(pcsource),1,out);
		}
	}
	printf("%d\n",Cnt);
	#ifndef USE_DISK
	rewind(out);
	fwrite(&NumSamples,sizeof(int),1,out);
	fwrite(&Cnt,sizeof(int),1,out);
	#endif

	fclose(in);
	fclose(out);

	return Cnt;
}


	
	
int convertAng(char infile[],int NumSamples){

/* Reads in data from file infile, converts it to x,y,z and outputs it to infile.bin.		*
 * Also prints information about the number of jackknife samples and number of data points.	*/

	int Cnt;
	char outfile[BUFFER_SIZE];
	double raRad,decRad,deg2rad,cdecRad;
	source sdata;
	pcsource pcdata;
	FILE *in,*out;
	#ifndef USE_DISK
	int i;
	#endif
	
	in = fopen(infile,"r");
	#ifdef USE_DISK
	sprintf(outfile,"%s.bin.temp",infile);
	#else
	sprintf(outfile,"%s0.bin",infile);
	#endif
	out = fopen(outfile,"w");
	Cnt=NumSamples;
	
	#ifndef USE_DISK
	fwrite(&Cnt,sizeof(int),1,out);
	for(i=0;i<=NumSamples;i++){
		fwrite(&Cnt,sizeof(int),1,out);
	}
	#endif
	deg2rad = M_PI/180.;
	Cnt = 0;
	while(!feof(in)){
		fscanf(in,"%lf %lf%*[^\n]",&sdata.ra,&sdata.dec);
		if(!feof(in)){
			Cnt++;
			raRad = deg2rad*sdata.ra;
			decRad = deg2rad*sdata.dec;
			cdecRad = cos(decRad);
			pcdata.x = cos(raRad)*cdecRad;
			pcdata.y = sin(raRad)*cdecRad;
			pcdata.z = sin(decRad);
			fwrite(&pcdata,sizeof(pcsource),1,out);
		}
	}
	printf("%d\n",Cnt);
	#ifndef USE_DISK
	rewind(out);
	fwrite(&NumSamples,sizeof(int),1,out);
	fwrite(&Cnt,sizeof(int),1,out);
	#endif

	fclose(in);
	fclose(out);

	return Cnt;
}

