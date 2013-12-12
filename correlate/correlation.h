#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "params.h"
#include "timing.h"
//#include <sys/types.h>
//#include <sys/stat.h>
//#include <fcntl.h>
//#include <unistd.h>

#ifdef USE_MPI
	#include "mpi.h"
#endif
#ifdef USE_OMP
	#include "omp.h"
#endif

/***********
    Macros
 ***********/

#ifdef DEBUGME
	#define DEBUGPRINT(string) fprintf(stdout,string);
#else
	#define DEBUGPRINT(string)
#endif

#define BINBUILD(bins) bins = (AngOrSpa==0) ? binBuildAngular(minDist,maxDist) : binBuildSpatial(minDist,maxDist)

#if defined (USE_MPI) || defined (USE_OMP)
	#define WORKNODES(tree,nodes) (AngOrSpa==0) ?												\
		angTreeWorkNodes((angTreeNode *)tree,WorkLevel,0,(angTreeNode **)nodes) :		\
		spaTreeWorkNodes((spaTreeNode *)tree,WorkLevel,0,(spaTreeNode **)nodes)
#else
	#define WORKNODES(tree,nodes)
#endif

/* setup AC to autocorr. data and CC to cross-corr. data using the appropriate routine */
#if defined (USE_MPI) && defined (USE_OMP)
	#define AC(srcs,tree1,tree2,bins) (AngOrSpa==0) ?												\
		ac_dist_shared(srcs,tree1,tree2,bins,NumWorkNodes,MyId) :		\
		ac_dist_shared(srcs,tree1,tree2,bins,NumWorkNodes,MyId)
	#define CC(srcs,targs,worknodes,tree1,tree2,bins) (AngOrSpa==0) ?									\
		cc_dist_shared(srcs,targs,worknodes,tree2,bins,NumWorkNodes,MyId) :	\
		cc_dist_shared(srcs,targs,worknodes,tree2,bins,NumWorkNodes,MyId)
#elif defined (USE_MPI)
	#define AC(sources,tree1,tree2,bins) (AngOrSpa==0) ?					 	\
		ac_dist(sources,tree1,tree2,bins,NumWorkNodes) :			\
		ac_dist(sources,tree1,tree2,bins,NumWorkNodes)
	#define CC(srcs,targs,worknodes,tree1,tree2,bins) (AngOrSpa==0) ?			 	\
		cc_dist(srcs,targs,worknodes,tree2,bins,NumWorkNodes) :		\
		cc_dist(srcs,targs,worknodes,tree2,bins,NumWorkNodes)
#elif defined (USE_OMP)
	#define AC(sources,tree1,tree2,bins) (AngOrSpa==0) ?						\
		ac_shared(sources,tree1,tree2,bins,NumWorkNodes,MyId) :		\
		ac_shared(sources,tree1,tree2,bins,NumWorkNodes,MyId)
	#define CC(srcs,targs,worknodes,tree1,tree2,bins) (AngOrSpa==0) ?				\
		cc_shared(srcs,targs,worknodes,tree2,bins,NumWorkNodes,MyId) :	\
		cc_shared(srcs,targs,worknodes,tree2,bins,NumWorkNodes,MyId)
#else
	#define AC(sources,tree1,tree2,bins) (AngOrSpa==0) ?						\
		ac_serial(sources,tree2,tree2,bins) :					\
		ac_serial(sources,tree2,tree2,bins)
	#define CC(srcs,targs,worknodes,tree1,tree2,bins) (AngOrSpa==0) ?				\
		cc_serial(srcs,targs,tree1,tree2,bins) :				\
		cc_serial(srcs,targs,tree1,tree2,bins)
#endif


/****************
 Global Variables
 ****************/

extern int NumBins,NumSamples,AngOrSpa;
#ifdef USE_MPI
extern int MyRank,NumProcs;
#endif
#ifdef USE_OMP
extern int NumThreads;
#endif


/****************************
 Global Structure Definitions
 ****************************/

typedef struct _pcsource{
	double x;
	double y;
	double z;
}pcsource;

typedef struct _bin{
	double limit;
	unsigned long long int *Cnt;
	double center;
}bin;

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

/*******************
 FUNCTION PROTOTYPES
********************/

/* input parameter function */
void getParams(char paramFile[],char fileList[],int *getDD,int *getDR,int *getRR,int *WorkLevel,double *minAngle,double *maxAngle);


/* correlation functions */
#ifdef USE_MPI
	#ifdef USE_OMP
		void ac_dist_shared(const pcsource data[],void *worknodes[],void *root,bin bins[],int NumWorkNodes,int MyId);
		void cc_dist_shared(const pcsource data[],const pcsource rand[],void *worknodes[],void *rootRandTree,bin bins[],int NumWorkNodes,int MyId);
	#else
		void ac_dist(const pcsource data[],void *worknodes[],void *root,bin bins[],int NumWorkNodes);
		void cc_dist(const pcsource data[],const pcsource rand[],void *worknodes[],void *rootRandTree,bin bins[],int NumWorkNodes);
	#endif
#else
	#ifdef USE_OMP
		void ac_shared(const pcsource data[],void *worknodes[],void *root,bin bins[],int NumWorkNodes,int MyId);
		void cc_shared(const pcsource data[],const pcsource rand[],void *worknodes[],void *rootRandTree,bin bins[],int NumWorkNodes,int MyId);
	#else
		void ac_serial(const pcsource data[],void *root1,void *root2,bin bins[]);
		void cc_serial(const pcsource dataA[],const pcsource dataB[],void *root1,void *root2,bin bins[]);
	#endif
#endif

/* pcsource functions */
void pcsourceGetMemory(pcsource **block1,pcsource **block2,char *fileList);
pcsource *pcsourceRead(char infile[],pcsource *block,int **Samples);

/* tree functions */
void *treeRead(char treeFile[], void *tblock);
void treeGetMemory(void **tblock1,void **tblock2,char *filelist);
#if defined (USE_MPI) || defined (USE_OMP)
int angTreeWorkNodes(angTreeNode *node,int WorkLevel,int CurNum,angTreeNode *worknodes[]);
int spaTreeWorkNodes(spaTreeNode *node,int WorkLevel,int CurNum,spaTreeNode *worknodes[]);
#endif

/* bin functions */
#ifdef USE_MPI
void binReduce(bin bins[]);
#endif
#ifdef USE_OMP
void binReduceOmp(bin bins[]);
#endif
bin *binBuildAngular(double minAngle,double maxAngle);
bin *binBuildSpatial(double minAngle,double maxAngle);
void binPrint(char filename[],bin bins[],int Samples1[],int Samples2[],int BinType);
void binFree(bin *bins);
void binClear(bin bins[]);

