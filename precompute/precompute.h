
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "params.h"

#ifdef USE_DISK
	#include <sys/types.h>
	#include <sys/mman.h>
	#include <unistd.h>
	#include <fcntl.h>
#endif

#define ERR_TOL 1.e-7

typedef struct _source{
	double ra;
	double dec;
}source;

typedef struct _zsource{
	double ra;
	double dec;
	double z;
}zsource;

typedef struct _pcsource{
	double x;
	double y;
	double z;
}pcsource;

typedef struct _spaTreeNode{
	int Sample;
	int Cnt;
	int Start;
	int End;
	double xmin,xmax,ymin,ymax,zmin,zmax;
	struct _spaTreeNode *lptr;
	struct _spaTreeNode *rptr;
}spaTreeNode;

typedef struct _angTreeNode{
	int Sample;
	int Cnt;
	int Start;
	int End;
	double x,y,z;
	double clowbound,c2lowbound,slowbound;
	struct _angTreeNode *lptr;
	struct _angTreeNode *rptr;
}angTreeNode;

typedef struct _splitlog{
	int dim;
	double val;
	int Cnt;
}splitlog;

typedef struct _filelist{
	int Start,End,Cnt;
	FILE *out;
	struct _filelist *next;
}filelist;

int convertAng(char infile[],int NumSamples);
int convertSpa(char infile[],int NumSamples);
int angTree(char filename[],splitlog dsplit[],int NumData,int NumSamples,int flag);
int spaTree(char filename[],splitlog dsplit[],int NumData,int NumSamples,int flag);
