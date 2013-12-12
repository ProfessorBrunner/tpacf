#include "precompute.h"

/* Useful macros */
#define SWAP(a,b) {     \
	Temp = (a);     \
	(a) = (b);      \
	(b) = Temp;     }
#define HUGEDBL 9999999999.

/* Some useful file scope variables */
static int isplit,CurSamp,CurFile;
static filelist *files,*hlist;

/******************************************************
 Contains routines to build angular/spatial trees and
 output results to file.
 ******************************************************/


/* Utility functions for both angular and spatial */

static double dividex(pcsource data[],int TreeInd[],int Start,int End,int *DivideInd){

	unsigned int i,j,Mid,k,Temp;
	double a;

	k = (Start+End)/2;
	*DivideInd = k;

	for( ; ; ){
		if(End <= Start+1){
			if(End == Start+1 && data[TreeInd[End]].x < data[TreeInd[Start]].x){
				SWAP(TreeInd[Start],TreeInd[End])
			}
			return data[TreeInd[k]].x;
		}else{
			Mid=(Start+End) >> 1;
			SWAP(TreeInd[Start+1],TreeInd[Mid])			
			if(data[TreeInd[Start]].x > data[TreeInd[End]].x){
				SWAP(TreeInd[Start],TreeInd[End])
			}
			if(data[TreeInd[Start+1]].x > data[TreeInd[End]].x){
				SWAP(TreeInd[Start+1],TreeInd[End])
			}
			if(data[TreeInd[Start]].x > data[TreeInd[Start+1]].x){
				SWAP(TreeInd[Start],TreeInd[Start+1])
			}
			i=Start+1;
			j=End;
			a=data[TreeInd[Start+1]].x;
			for( ; ; ){
				do{
					i++;
				}while (data[TreeInd[i]].x < a);
				do{
					j--;
				}while (data[TreeInd[j]].x > a);
				if (j < i) break;
				SWAP(TreeInd[i],TreeInd[j])
			}
			SWAP(TreeInd[Start+1],TreeInd[j])
			if(j >= k){
				End=j-1;
			}
			if(j <= k){
				Start=i;
			}
		}
	}
}

static double dividey(pcsource data[],int TreeInd[],int Start,int End,int *DivideInd){

	unsigned int i,j,Mid,k,Temp;
	double a;

	k = (Start+End)/2;
	*DivideInd = k;

	for( ; ; ){
		if(End <= Start+1){
			if(End == Start+1 && data[TreeInd[End]].y < data[TreeInd[Start]].y){
				SWAP(TreeInd[Start],TreeInd[End])
			}
			return data[TreeInd[k]].y;
		}else{
			Mid=(Start+End) >> 1;
			SWAP(TreeInd[Start+1],TreeInd[Mid])			
			if(data[TreeInd[Start]].y > data[TreeInd[End]].y){
				SWAP(TreeInd[Start],TreeInd[End])
			}
			if(data[TreeInd[Start+1]].y > data[TreeInd[End]].y){
				SWAP(TreeInd[Start+1],TreeInd[End])
			}
			if(data[TreeInd[Start]].y > data[TreeInd[Start+1]].y){
				SWAP(TreeInd[Start],TreeInd[Start+1])
			}
			i=Start+1;
			j=End;
			a=data[TreeInd[Start+1]].y;
			for( ; ; ){
				do{
					i++;
				}while (data[TreeInd[i]].y < a);
				do{
					j--;
				}while (data[TreeInd[j]].y > a);
				if(j < i) break;
				SWAP(TreeInd[i],TreeInd[j])
			}
			SWAP(TreeInd[Start+1],TreeInd[j])
			if(j >= k){
				End=j-1;
			}
			if(j <= k){
				Start=i;
			}
		}
	}
}

static double dividez(pcsource data[],int TreeInd[],int Start,int End,int *DivideInd){

	unsigned int i,j,Mid,k,Temp;
	double a;
	
	k = (Start+End)/2;
	*DivideInd = k;

	for( ; ; ){
		if(End <= Start+1){
			if(End == Start+1 && data[TreeInd[End]].z < data[TreeInd[Start]].z){
				SWAP(TreeInd[Start],TreeInd[End])
			}
			return data[TreeInd[k]].z;
		}else{
			Mid=(Start+End) >> 1;
			SWAP(TreeInd[Start+1],TreeInd[Mid])			
			if(data[TreeInd[Start]].z > data[TreeInd[End]].z){
				SWAP(TreeInd[Start],TreeInd[End])
			}
			if(data[TreeInd[Start+1]].z > data[TreeInd[End]].z){
				SWAP(TreeInd[Start+1],TreeInd[End])
			}
			if(data[TreeInd[Start]].z > data[TreeInd[Start+1]].z){
				SWAP(TreeInd[Start],TreeInd[Start+1])
			}
			i=Start+1;
			j=End;
			a=data[TreeInd[Start+1]].z;
			for( ; ; ){
				do{
					i++;
				}while (data[TreeInd[i]].z < a);
				do{
					j--;
				}while (data[TreeInd[j]].z > a);
				if (j < i) break;
				SWAP(TreeInd[i],TreeInd[j])
			}
			SWAP(TreeInd[Start+1],TreeInd[j])
			if(j >= k){
				End=j-1;
			}
			if(j <= k){
				Start=i;
			}
		}
	}
}

static void dividexx(pcsource data[],int TreeInd[],int Start,int End,double split,int *DivideInd){

	int i,j,Temp;
	
	i = Start;
	j = End;

	while(i<j){

		while(data[TreeInd[i]].x <= split && i < j){
			i++;
		}
	
		while(data[TreeInd[j]].x > split && j > i){
			j--;
		}

		if(i != j){
			SWAP(TreeInd[i],TreeInd[j])
		}else{
			if(data[TreeInd[i]].x > split){
				i--;
			}
			break;
		}
	}

	*DivideInd = i;

	return;
}

static void divideyy(pcsource data[],int TreeInd[],int Start,int End,double split,int *DivideInd){

	int i,j,Temp;

	i = Start;
	j = End;

	while(i<j){

		while(data[TreeInd[i]].y <= split && i < j){
			i++;
		}
	
		while(data[TreeInd[j]].y > split && j > i){
			j--;
		}

		if(i != j){
			SWAP(TreeInd[i],TreeInd[j])
		}else{
			if(data[TreeInd[i]].y > split){
				i--;
			}
			break;
		}
	}

	*DivideInd = i;

	return;
}

static void dividezz(pcsource data[],int TreeInd[],int Start,int End,double split,int *DivideInd){

	int i,j,Temp;

	i = Start;
	j = End;

	while(i<j){

		while(data[TreeInd[i]].z <= split && i < j){
			i++;
		}
	
		while(data[TreeInd[j]].z > split && j > i){
			j--;
		}

		if(i != j){
			SWAP(TreeInd[i],TreeInd[j])
		}else{
			if(data[TreeInd[i]].z > split){
				i--;
			}
			break;
		}
	}

	*DivideInd = i;

	return;
}



/***************** Spatial Tree functions ******************/

static void spaTreePrint(spaTreeNode *node, int Leaf, FILE *out){

	if(Leaf == 1){
		node->lptr = (spaTreeNode *)1;	/* Encode that the node has children */
	}
	fwrite(node,sizeof(spaTreeNode),1,out);
	
	return;
}

static void spaTreeBuild(spaTreeNode *node,pcsource data[],int TreeInd[],int LogSamples,int Depth,splitlog dsplit[],char dataname[],FILE *out){

/* Builds tree for spatial data */

	char filename[BUFFER_SIZE];
	int i,Dim,DivideInd,Cnt,Sample,Start,End;
	double split,xsize,ysize,zsize;
	double smallx,smally,smallz,bigx,bigy,bigz;
	spaTreeNode lNewNode,rNewNode;
	filelist *tfiles;
	typedef struct _stack{
		spaTreeNode node;
		int Depth;
		FILE *myfile;
		struct _stack *next;
	}stack;
	stack *nstack,*nnstack;


	nstack = malloc(sizeof(stack));
	if(nstack==NULL){
		fprintf(stderr,"Failed to allocate stack\n");
		return;
	}
	memcpy(&(nstack->node),node,sizeof(spaTreeNode));
	nstack->Depth = Depth;
	nstack->myfile = out;
	nstack->next = NULL;

	while(nstack != NULL){

		Cnt = (nstack->node).Cnt;
		Sample = (nstack->node).Sample;
		Start = (nstack->node).Start;
		End = (nstack->node).End;
		Depth = nstack->Depth;

		if((nstack->node).Cnt <= MAXNUM_PER_FILE){
			if(nstack->myfile == NULL){
				if(Depth == 0){
					sprintf(filename,"%s0.tree",dataname);
					if((hlist = malloc(sizeof(filelist))) == NULL){
						fprintf(stderr,"Failed to allocate tfiles\n");
						return;
					}
					files = hlist;
					files->Start = 0;
					files->Cnt = 0;
					CurFile++;
				}else{
					if((tfiles = malloc(sizeof(filelist))) == NULL){
						fprintf(stderr,"Failed to allocate tfiles\n");
						return;
					}
					tfiles->Start = Start;
					tfiles->End = End;
					tfiles->Cnt = 0;
					tfiles->next = NULL;
					if(CurFile == 0){
						hlist = tfiles;
					}else{
						files->next = tfiles;
					}
					files = tfiles;
					sprintf(filename,"%s%d.tree",dataname,CurFile);
					CurFile++;
				}
				nstack->myfile = fopen(filename,"w");
				files->out = nstack->myfile;
				fwrite(&CurFile,sizeof(int),1,nstack->myfile);
			}
			if((nstack->node).Cnt > MAX_SPATIAL_NODE_CNT){
				(nstack->node).Start -= files->Start;
				(nstack->node).End -= files->Start;
				spaTreePrint(&(nstack->node),1,nstack->myfile);
				files->Cnt++;
			}else{
				(nstack->node).Start -= files->Start;
				(nstack->node).End -= files->Start;
				spaTreePrint(&(nstack->node),0,nstack->myfile);
				files->Cnt++;
				nnstack = nstack;
				nstack = nstack->next;
				free(nnstack);
				continue;
			}
		}
		
		xsize = (nstack->node).xmax - (nstack->node).xmin;
		ysize = (nstack->node).ymax - (nstack->node).ymin;
		zsize = (nstack->node).zmax - (nstack->node).zmin;
		if(xsize > ysize){
			if(xsize > zsize){
				Dim = 0;
				split = dividex(data,TreeInd,Start,End,&DivideInd);
			}else{
				Dim = 2;
				split = dividez(data,TreeInd,Start,End,&DivideInd);
			}
		}else{
			if(ysize > zsize){
				Dim = 1;
				split = dividey(data,TreeInd,Start,End,&DivideInd);
			}else{
				Dim = 2;
				split = dividez(data,TreeInd,Start,End,&DivideInd);
			}
		}

		if(Depth < LogSamples){
			dsplit[isplit].dim = Dim;
			dsplit[isplit].val = split;
			isplit++;
		}
		
		lNewNode.Sample = Sample;
		lNewNode.Cnt = DivideInd-Start+1;
		lNewNode.Start = Start;
		lNewNode.End = DivideInd;
		lNewNode.lptr = NULL;
		lNewNode.rptr = NULL;
		
		rNewNode.Sample = Sample;
		rNewNode.Cnt = Cnt - lNewNode.Cnt;
		rNewNode.Start = DivideInd+1;
		rNewNode.End = End;
		rNewNode.lptr = NULL;
		rNewNode.rptr = NULL;

		if(Depth == LogSamples-1){
			dsplit[CurSamp-1].Cnt = lNewNode.Cnt;
			dsplit[CurSamp].Cnt = rNewNode.Cnt;
			lNewNode.Sample = CurSamp;
			CurSamp++;
			rNewNode.Sample = CurSamp;
			CurSamp++;
		}
		
		smallx = smally = smallz = HUGEDBL;
		bigx = bigy = bigz = -HUGEDBL;

		for(i=lNewNode.Start;i<=lNewNode.End;i++){
			if(data[TreeInd[i]].x > bigx) bigx = data[TreeInd[i]].x;
			if(data[TreeInd[i]].x < smallx) smallx = data[TreeInd[i]].x;
			if(data[TreeInd[i]].y > bigy) bigy = data[TreeInd[i]].y;
			if(data[TreeInd[i]].y < smally) smally = data[TreeInd[i]].y;
			if(data[TreeInd[i]].z > bigz) bigz = data[TreeInd[i]].z;
			if(data[TreeInd[i]].z < smallz) smallz = data[TreeInd[i]].z;
		}

		lNewNode.xmin = smallx;
		lNewNode.xmax = bigx;
		lNewNode.ymin = smally;
		lNewNode.ymax = bigy;
		lNewNode.zmin = smallz;
		lNewNode.zmax = bigz;
	
		smallx = smally = smallz = HUGEDBL;
		bigx = bigy = bigz = -HUGEDBL;

		for(i=rNewNode.Start;i<=rNewNode.End;i++){
			if(data[TreeInd[i]].x > bigx) bigx = data[TreeInd[i]].x;
			if(data[TreeInd[i]].x < smallx) smallx = data[TreeInd[i]].x;
			if(data[TreeInd[i]].y > bigy) bigy = data[TreeInd[i]].y;
			if(data[TreeInd[i]].y < smally) smally = data[TreeInd[i]].y;
			if(data[TreeInd[i]].z > bigz) bigz = data[TreeInd[i]].z;
			if(data[TreeInd[i]].z < smallz) smallz = data[TreeInd[i]].z;
		}

		rNewNode.xmin = smallx;
		rNewNode.xmax = bigx;
		rNewNode.ymin = smally;
		rNewNode.ymax = bigy;
		rNewNode.zmin = smallz;
		rNewNode.zmax = bigz;

		memcpy(&(nstack->node),&rNewNode,sizeof(spaTreeNode));
		nstack->Depth = Depth+1;
		nnstack = malloc(sizeof(stack));
		if(nnstack == NULL){
			fprintf(stderr,"Failed to allocate stack space in buildtree\n");
			return;
		}
		memcpy(&(nnstack->node),&lNewNode,sizeof(spaTreeNode));
		nnstack->Depth = Depth+1;
		nnstack->myfile = nstack->myfile;
		nnstack->next = nstack;
		nstack = nnstack;
	}

	return; 
}

static void spaTreeBuildRan(spaTreeNode *node,pcsource rand[],int TreeInd[],int LogSamples,splitlog dsplit[],int Depth,char dataname[],FILE *out){

/* Builds tree for data, based on the tree built for a previous data set.  Designed to ensure that all jackknife	*
 * samples occupy the same regions across different data sets in the same calculation.					*/

 	char filename[BUFFER_SIZE];
	int i,DivideInd;
	double smallx,smally,smallz,bigx,bigy,bigz;
	spaTreeNode *rNewNode,*lNewNode;
	spaTreeNode ncpy;
	filelist *tfiles;

	if(Depth == LogSamples){
		node->Sample = CurSamp;
		dsplit[CurSamp-1].Cnt = node->Cnt;
		CurSamp++;
		spaTreeBuild(node,rand,TreeInd,LogSamples,Depth,dsplit,dataname,out);
		return;
	}

	if(node->Cnt <= MAXNUM_PER_FILE){
		if(out != NULL){
			memcpy(&ncpy,node,sizeof(spaTreeNode));
			ncpy.Start -= files->Start;
			ncpy.End -= files->Start;
			spaTreePrint(&ncpy,1,out);
			files->Cnt++;
		}else{
			if(Depth == 0){
				sprintf(filename,"%s0.tree",dataname);
				if((hlist = malloc(sizeof(filelist))) == NULL){
					fprintf(stderr,"Failed to allocate tfiles\n");
					return;
				}
				files = hlist;
				files->Start = 0;
				files->Cnt = 0;
				CurFile++;
			}else{
				if((tfiles = malloc(sizeof(filelist))) == NULL){
					fprintf(stderr,"Failed to allocate tfiles\n");
					return;
				}
				tfiles->Start = node->Start;
				tfiles->End = node->End;
				tfiles->Cnt = 0;
				tfiles->next = NULL;
				if(CurFile == 0){
					hlist = tfiles;
				}else{
					files->next = tfiles;
				}
				files = tfiles;
				sprintf(filename,"%s%d.tree",dataname,CurFile);
				CurFile++;
			}
			out = fopen(filename,"w");
			fwrite(&CurFile,sizeof(int),1,out);
			files->out = out;
			memcpy(&ncpy,node,sizeof(spaTreeNode));
			ncpy.Start -= files->Start;
			ncpy.End -= files->Start;
			spaTreePrint(&ncpy,1,out);
			files->Cnt++;
		}
	}
	
	lNewNode = malloc(sizeof(spaTreeNode));
	rNewNode = malloc(sizeof(spaTreeNode));
	if(lNewNode == NULL || rNewNode == NULL){
		printf("Failed to allocate memory in buildtreeran...\n");
		return;
	}

	node->lptr = lNewNode;
	node->rptr = rNewNode;

	if(dsplit[isplit].dim == 0){
		dividexx(rand,TreeInd,node->Start,node->End,dsplit[isplit].val,&DivideInd);
	}else if(dsplit[isplit].dim == 1){
		divideyy(rand,TreeInd,node->Start,node->End,dsplit[isplit].val,&DivideInd);
	}else{
		dividezz(rand,TreeInd,node->Start,node->End,dsplit[isplit].val,&DivideInd);
	}
	isplit++;

	lNewNode->Sample = node->Sample;
	lNewNode->Cnt = DivideInd - node->Start + 1;
	lNewNode->Start = node->Start;
	lNewNode->End = DivideInd;
	lNewNode->lptr = NULL;
	lNewNode->rptr = NULL;
	rNewNode->Sample = node->Sample;
	rNewNode->Cnt = node->Cnt - lNewNode->Cnt;
	rNewNode->Start = DivideInd+1;
	rNewNode->End = node->End;
	rNewNode->lptr = NULL;
	rNewNode->rptr = NULL;

	smallx = smally = smallz = HUGEDBL;
	bigx = bigy = bigz = -HUGEDBL;

	for(i=lNewNode->Start;i<=lNewNode->End;i++){
		if(rand[TreeInd[i]].x > bigx) bigx = rand[TreeInd[i]].x;
		if(rand[TreeInd[i]].x < smallx) smallx = rand[TreeInd[i]].x;
		if(rand[TreeInd[i]].y > bigy) bigy = rand[TreeInd[i]].y;
		if(rand[TreeInd[i]].y < smally) smally = rand[TreeInd[i]].y;
		if(rand[TreeInd[i]].z > bigz) bigz = rand[TreeInd[i]].z;
		if(rand[TreeInd[i]].z < smallz) smallz = rand[TreeInd[i]].z;
	}

	lNewNode->xmin = smallx;
	lNewNode->xmax = bigx;
	lNewNode->ymin = smally;
	lNewNode->ymax = bigy;
	lNewNode->zmin = smallz;
	lNewNode->zmax = bigz;

	smallx = smally = smallz = HUGEDBL;
	bigx = bigy = bigz = -HUGEDBL;

	for(i=rNewNode->Start;i<=rNewNode->End;i++){
		if(rand[TreeInd[i]].x > bigx) bigx = rand[TreeInd[i]].x;
		if(rand[TreeInd[i]].x < smallx) smallx = rand[TreeInd[i]].x;
		if(rand[TreeInd[i]].y > bigy) bigy = rand[TreeInd[i]].y;
		if(rand[TreeInd[i]].y < smally) smally = rand[TreeInd[i]].y;
		if(rand[TreeInd[i]].z > bigz) bigz = rand[TreeInd[i]].z;
		if(rand[TreeInd[i]].z < smallz) smallz = rand[TreeInd[i]].z;
	}
	
	rNewNode->xmin = smallx;
	rNewNode->xmax = bigx;
	rNewNode->ymin = smally;
	rNewNode->ymax = bigy;
	rNewNode->zmin = smallz;
	rNewNode->zmax = bigz;

	spaTreeBuildRan(lNewNode,rand,TreeInd,LogSamples,dsplit,Depth+1,dataname,out);
	spaTreeBuildRan(rNewNode,rand,TreeInd,LogSamples,dsplit,Depth+1,dataname,out);

	return;
}


int spaTree(char filename[],splitlog dsplit[],int NumData,int NumSamples,int flag){

/* Coordinates tree building process, printing the results, and freeing memory allocated for a tree.	*
 * Also reorganizes data based on the tree and prints sample counts in the ".bin" data file.		*/

	int i,j,*TreeInd,LogSamples;
	char dfile[BUFFER_SIZE];
	double smallx,bigx,smally,bigy,smallz,bigz;
	spaTreeNode *root;
	pcsource *data;
	FILE *inout;
	#ifdef USE_DISK
		int pin;
	#endif

	if(NumSamples != 0){
		LogSamples = (int)(log((double)NumSamples)/M_LN2+0.5);
	}else{
		LogSamples = 0;
	}

	root = malloc(sizeof(spaTreeNode));
	if(root == NULL){
		fprintf(stderr,"Failed to allocate root node\n");
		return 0;
	}

	#ifdef USE_DISK
		sprintf(dfile,"%s.bin.temp",filename);
		pin = open(dfile,O_RDONLY);
		data = mmap(NULL,NumData*sizeof(pcsource),PROT_READ,MAP_SHARED,pin,0);
		if((int)data == -1){
			fprintf(stderr,"Failed to map file to memory\n");
			return 0;
		}
		close(pin);
	#else
		sprintf(dfile,"%s0.bin",filename);
		inout = fopen(dfile,"r");
		fread(&i,sizeof(int),1,inout);
		fread(&NumData,sizeof(int),1,inout);
		for(i=0;i<NumSamples;i++){
			fread(&dsplit[i].Cnt,sizeof(int),1,inout);
		}

		data = malloc(NumData*sizeof(pcsource));
		if(data == NULL){
			fprintf(stderr,"Failled to allocate data in tree\n");
			return 0;
		}
		fread(data,sizeof(pcsource),NumData,inout);
		fclose(inout);
		remove(dfile);
	#endif
	
	TreeInd = malloc(NumData*sizeof(int));
	if(TreeInd == NULL){
		fprintf(stderr,"Failed to allocate TreeInd\n");
		return 0;
	}

	root->Sample = -1;
	root->Cnt = NumData;
	root->Start = 0;
	root->End = NumData-1;
	root->lptr = NULL;
	root->rptr = NULL;

	smallx = smally = smallz = HUGEDBL;
	bigx = bigy = bigz = -HUGEDBL;

	for(i=0;i<NumData;i++){
		TreeInd[i] = i;
		if(data[i].x > bigx) bigx = data[i].x;
		if(data[i].x < smallx) smallx = data[i].x;
		if(data[i].y > bigy) bigy = data[i].y;
		if(data[i].y < smally) smally = data[i].y;
		if(data[i].z > bigz) bigz = data[i].z;
		if(data[i].z < smallz) smallz = data[i].z;
	}

	root->xmin = smallx;
	root->xmax = bigx;
	root->ymin = smally;
	root->ymax = bigy;
	root->zmin = smallz;
	root->zmax = bigz;

	isplit = 0;
	CurSamp = 1;
	CurFile = 0;
	if(flag == 0){
		spaTreeBuild(root,data,TreeInd,LogSamples,0,dsplit,filename,NULL);
	}else{
		spaTreeBuildRan(root,data,TreeInd,LogSamples,dsplit,0,filename,NULL);
	}
	
	if(CurFile == 1){
		rewind(hlist->out);
		fwrite(&(hlist->Cnt),sizeof(int),1,hlist->out);
		fclose(hlist->out);
		if(root != NULL) free(root);
		sprintf(dfile,"%s0.bin",filename);
		inout = fopen(dfile,"w");
		fwrite(&NumSamples,sizeof(int),1,inout);
		fwrite(&NumData,sizeof(int),1,inout);
		for(i=0;i<NumSamples;i++){
			fwrite(&dsplit[i].Cnt,sizeof(int),1,inout);
		}
		for(i=0;i<NumData;i++){
			fwrite(&(data[TreeInd[i]].x),sizeof(pcsource),1,inout);
		}
		fclose(inout);
	}else{
		if(root != NULL) free(root);
		for(i=0;i<CurFile;i++){
			rewind(hlist->out);
			fwrite(&(hlist->Cnt),sizeof(int),1,hlist->out);
			fclose(hlist->out);
			sprintf(dfile,"%s%d.bin",filename,i);
			inout = fopen(dfile,"w");
			NumData = hlist->End - hlist->Start + 1;
			fwrite(&NumSamples,sizeof(int),1,inout);
			fwrite(&NumData,sizeof(int),1,inout);
			for(j=0;j<NumSamples;j++){
				fwrite(&dsplit[j].Cnt,sizeof(int),1,inout);
			}
			for(j=hlist->Start;j<=hlist->End;j++){
				fwrite(&(data[TreeInd[j]].x),sizeof(pcsource),1,inout);
			}
			fclose(inout);
			files = hlist->next;
			free(hlist);
			hlist = files;
		}
	}	

	free(TreeInd);
	
	#ifdef USE_DISK
	munmap(data,NumData*sizeof(pcsource));
	sprintf(dfile,"%s.bin.temp",filename);
	remove(dfile);
	#else
	fflush(stdout);
	free(data);
	#endif

	return CurFile;
}



/*********** Angular tree functions ********************/

static void angTreePrint(angTreeNode *node, int Leaf, FILE *out){

	if(Leaf == 1){
		node->lptr = (angTreeNode *)1;		/* Encode that the node has children */
	}
	fwrite(node,sizeof(angTreeNode),1,out);
	
	return;
}

static void angTreeBuild(angTreeNode *node,double x[],double y[],double z[],pcsource data[],int TreeInd[],int Depth,int LogSamples,splitlog dsplit[],char dataname[],FILE *out){

/* Builds tree for angular data */

	char filename[BUFFER_SIZE];
	int i,Dim,DivideInd,Cnt,Sample,Start,End;
	double split,xsize,ysize,zsize,norm;
	double smallx,smally,smallz,bigx,bigy,bigz,avgx,avgy,avgz;
	double mind,mindsq,d;
	angTreeNode lNewNode,rNewNode;
	filelist *tfiles;
	typedef struct _stack{
		angTreeNode node;
		double x[2];
		double y[2];
		double z[2];
		int Depth;
		FILE *myfile;
		struct _stack *next;
	}stack;
	stack *nstack,*nnstack;
	
	nstack = malloc(sizeof(stack));
	if(nstack == NULL){
		printf("Failed to allocate stack space in buildtree\n");
		return;
	}
	memcpy(&(nstack->node),node,sizeof(angTreeNode));
	nstack->x[0] = x[0];
	nstack->x[1] = x[1];
	nstack->y[0] = y[0];
	nstack->y[1] = y[1];
	nstack->z[0] = z[0];
	nstack->z[1] = z[1];
	nstack->Depth = Depth;
	nstack->myfile = out;
	nstack->next = NULL;

	while(nstack != NULL){

		Cnt = (nstack->node).Cnt;
		Sample = (nstack->node).Sample;
		Start = (nstack->node).Start;
		End = (nstack->node).End;
		x[0] = nstack->x[0];
		x[1] = nstack->x[1];
		y[0] = nstack->y[0];
		y[1] = nstack->y[1];
		z[0] = nstack->z[0];
		z[1] = nstack->z[1];
		Depth = nstack->Depth;

		if((nstack->node).Cnt <= MAXNUM_PER_FILE){
			if(nstack->myfile == NULL){
				if(Depth == 0){
					sprintf(filename,"%s0.tree",dataname);
					if((hlist = malloc(sizeof(filelist))) == NULL){
						fprintf(stderr,"Failed to allocate tfiles\n");
						return;
					}
					files = hlist;
					files->Start = 0;
					files->Cnt = 0;
					CurFile++;
				}else{
					if((tfiles = malloc(sizeof(filelist))) == NULL){
						fprintf(stderr,"Failed to allocate tfiles\n");
						return;
					}
					tfiles->Start = Start;
					tfiles->End = End;
					tfiles->Cnt = 0;
					tfiles->next = NULL;
					if(CurFile == 0){
						hlist = tfiles;
					}else{
						files->next = tfiles;
					}
					files = tfiles;
					sprintf(filename,"%s%d.tree",dataname,CurFile);
					CurFile++;
				}
				nstack->myfile = fopen(filename,"w");
				files->out = nstack->myfile;
				fwrite(&CurFile,sizeof(int),1,nstack->myfile);
			}
			if((nstack->node).Cnt > MAX_ANGULAR_NODE_CNT){
				(nstack->node).Start -= files->Start;
				(nstack->node).End -= files->Start;
				angTreePrint(&(nstack->node),1,nstack->myfile);
				files->Cnt++;
			}else{
				(nstack->node).Start -= files->Start;
				(nstack->node).End -= files->Start;
				angTreePrint(&(nstack->node),0,nstack->myfile);
				files->Cnt++;
				nnstack = nstack;
				nstack = nstack->next;
				free(nnstack);
				continue;
			}
		}

		xsize = x[1] - x[0];
		ysize = y[1] - y[0];
		zsize = z[1] - z[0];
		if(xsize > ysize){
			if(xsize > zsize){
				Dim = 0;
				split = dividex(data,TreeInd,Start,End,&DivideInd);
			}else{
				Dim = 2;
				split = dividez(data,TreeInd,Start,End,&DivideInd);
			}
		}else{
			if(ysize > zsize){
				Dim = 1;
				split = dividey(data,TreeInd,Start,End,&DivideInd);
			}else{
				Dim = 2;
				split = dividez(data,TreeInd,Start,End,&DivideInd);
			}
		}

		if(Depth < LogSamples){
			dsplit[isplit].dim = Dim;
			dsplit[isplit].val = split;
			isplit++;
		}
		
		lNewNode.Sample = Sample;
		lNewNode.Cnt = DivideInd-Start+1;
		lNewNode.Start = Start;
		lNewNode.End = DivideInd;
		lNewNode.lptr = NULL;
		lNewNode.rptr = NULL;

		rNewNode.Sample = Sample;
		rNewNode.Cnt = Cnt - lNewNode.Cnt;
		rNewNode.Start = DivideInd+1;
		rNewNode.End = End;
		rNewNode.lptr = NULL;
		rNewNode.rptr = NULL;
		
		if(Depth == LogSamples-1){
			dsplit[CurSamp-1].Cnt = lNewNode.Cnt;
			dsplit[CurSamp].Cnt = rNewNode.Cnt;
			lNewNode.Sample = CurSamp;
			CurSamp++;
			rNewNode.Sample = CurSamp;
			CurSamp++;
		}
		
		smallx = smally = smallz = HUGEDBL;
		bigx = bigy = bigz = -HUGEDBL;

		for(i=lNewNode.Start;i<=lNewNode.End;i++){
			if(data[TreeInd[i]].x > bigx) bigx = data[TreeInd[i]].x;
			if(data[TreeInd[i]].x < smallx) smallx = data[TreeInd[i]].x;
			if(data[TreeInd[i]].y > bigy) bigy = data[TreeInd[i]].y;
			if(data[TreeInd[i]].y < smally) smally = data[TreeInd[i]].y;
			if(data[TreeInd[i]].z > bigz) bigz = data[TreeInd[i]].z;
			if(data[TreeInd[i]].z < smallz) smallz = data[TreeInd[i]].z;
		}

		x[0] = smallx;
		x[1] = bigx;
		y[0] = smally;
		y[1] = bigy;
		z[0] = smallz;
		z[1] = bigz;

		avgx = (smallx+bigx)/2.;
		avgy = (smally+bigy)/2.;
		avgz = (smallz+bigz)/2.;
		norm = sqrt(avgx*avgx + avgy*avgy + avgz*avgz);
		lNewNode.x = avgx/norm;
		lNewNode.y = avgy/norm;
		lNewNode.z = avgz/norm;

		mind = HUGEDBL;
		for(i=lNewNode.Start;i<=lNewNode.End;i++){
			d = lNewNode.x*data[TreeInd[i]].x + lNewNode.y*data[TreeInd[i]].y + lNewNode.z*data[TreeInd[i]].z;
			if(d < mind) mind = d;
		}
		lNewNode.clowbound = mind;
		mindsq = mind*mind;
		lNewNode.c2lowbound = 2.*mindsq-1.;
		lNewNode.slowbound = sqrt(1-mindsq);


		smallx = smally = smallz = HUGEDBL;
		bigx = bigy = bigz = -HUGEDBL;

		for(i=rNewNode.Start;i<=rNewNode.End;i++){
			if(data[TreeInd[i]].x > bigx) bigx = data[TreeInd[i]].x;
			if(data[TreeInd[i]].x < smallx) smallx = data[TreeInd[i]].x;
			if(data[TreeInd[i]].y > bigy) bigy = data[TreeInd[i]].y;
			if(data[TreeInd[i]].y < smally) smally = data[TreeInd[i]].y;
			if(data[TreeInd[i]].z > bigz) bigz = data[TreeInd[i]].z;
			if(data[TreeInd[i]].z < smallz) smallz = data[TreeInd[i]].z;
		}

		avgx = (smallx+bigx)/2.;
		avgy = (smally+bigy)/2.;
		avgz = (smallz+bigz)/2.;
		norm = sqrt(avgx*avgx + avgy*avgy + avgz*avgz);
		rNewNode.x = avgx/norm;
		rNewNode.y = avgy/norm;
		rNewNode.z = avgz/norm;

		mind = HUGEDBL;
		for(i=rNewNode.Start;i<=rNewNode.End;i++){
			d = rNewNode.x*data[TreeInd[i]].x + rNewNode.y*data[TreeInd[i]].y + rNewNode.z*data[TreeInd[i]].z;
			if(d < mind) mind = d;
		}
		rNewNode.clowbound = mind;
		mindsq = mind*mind;
		rNewNode.c2lowbound = 2.*mindsq-1.;
		rNewNode.slowbound = sqrt(1-mindsq);

		memcpy(&(nstack->node),&rNewNode,sizeof(angTreeNode));
		nstack->x[0] = smallx;
		nstack->x[1] = bigx;
		nstack->y[0] = smally;
		nstack->y[1] = bigy;
		nstack->z[0] = smallz;
		nstack->z[1] = bigz;
		nstack->Depth = Depth+1;
		nnstack = malloc(sizeof(stack));
		if(nnstack == NULL){
			fprintf(stderr,"Failed to allocate stack space in buildtree\n");
			return;
		}
		memcpy(&(nnstack->node),&lNewNode,sizeof(angTreeNode));
		nnstack->x[0] = x[0];
		nnstack->x[1] = x[1];
		nnstack->y[0] = y[0];
		nnstack->y[1] = y[1];
		nnstack->z[0] = z[0];
		nnstack->z[1] = z[1];
		nnstack->Depth = Depth+1;
		nnstack->myfile = nstack->myfile;
		nnstack->next = nstack;
		nstack = nnstack;
	}

	return; 
}

static void angTreeBuildRan(angTreeNode *node,double x[],double y[],double z[],pcsource rand[],int TreeInd[],int LogSamples,splitlog dsplit[],int Depth,char dataname[],FILE *out){

/* Builds tree for data, based on the tree built for a previous data set.  Designed to ensure that all jackknife	*
 * samples occupy the same regions across different data sets in the same calculation.					*/

 	char filename[BUFFER_SIZE];
	int i,DivideInd;
	double smallx,smally,smallz,bigx,bigy,bigz,avgx,avgy,avgz,norm,d,mind,mindsq;
	angTreeNode *rNewNode,*lNewNode;
	angTreeNode ncpy;
	filelist *tfiles;

	if(Depth == LogSamples){
		node->Sample = CurSamp;
		dsplit[CurSamp-1].Cnt = node->Cnt;
		CurSamp++;
		angTreeBuild(node,x,y,z,rand,TreeInd,Depth,LogSamples,dsplit,dataname,out);
		return;
	}
	
	if(node->Cnt <= MAXNUM_PER_FILE){
		if(out != NULL){
			memcpy(&ncpy,node,sizeof(angTreeNode));
			ncpy.Start -= files->Start;
			ncpy.End -= files->Start;
			angTreePrint(&ncpy,1,out);
			files->Cnt++;
		}else{
			if(Depth == 0){
				sprintf(filename,"%s0.tree",dataname);
				if((hlist = malloc(sizeof(filelist))) == NULL){
					fprintf(stderr,"Failed to allocate tfiles\n");
					return;
				}
				files = hlist;
				files->Start = 0;
				files->Cnt = 0;
				CurFile++;
			}else{
				if((tfiles = malloc(sizeof(filelist))) == NULL){
					fprintf(stderr,"Failed to allocate tfiles\n");
					return;
				}
				tfiles->Start = node->Start;
				tfiles->End = node->End;
				tfiles->Cnt = 0;
				tfiles->next = NULL;
				if(CurFile == 0){
					hlist = tfiles;
				}else{
					files->next = tfiles;
				}
				files = tfiles;
				sprintf(filename,"%s%d.tree",dataname,CurFile);
				CurFile++;
			}
			out = fopen(filename,"w");
			fwrite(&CurFile,sizeof(int),1,out);
			files->out = out;
			memcpy(&ncpy,node,sizeof(angTreeNode));
			ncpy.Start -= files->Start;
			ncpy.End -= files->Start;
			angTreePrint(&ncpy,1,out);
			files->Cnt++;
		}
	}

	lNewNode = malloc(sizeof(angTreeNode));
	rNewNode = malloc(sizeof(angTreeNode));
	if(lNewNode == NULL || rNewNode == NULL){
		printf("Failed to allocate memory in buildtreeran...\n");
		return;
	}

	node->lptr = lNewNode;
	node->rptr = rNewNode;

	if(dsplit[isplit].dim == 0){
		dividexx(rand,TreeInd,node->Start,node->End,dsplit[isplit].val,&DivideInd);
	}else if(dsplit[isplit].dim == 1){
		divideyy(rand,TreeInd,node->Start,node->End,dsplit[isplit].val,&DivideInd);
	}else{
		dividezz(rand,TreeInd,node->Start,node->End,dsplit[isplit].val,&DivideInd);
	}
	isplit++;

	lNewNode->Sample = node->Sample;
	lNewNode->Cnt = DivideInd - node->Start + 1;
	lNewNode->Start = node->Start;
	lNewNode->End = DivideInd;
	lNewNode->lptr = NULL;
	lNewNode->rptr = NULL;
	rNewNode->Sample = node->Sample;
	rNewNode->Cnt = node->Cnt - lNewNode->Cnt;
	rNewNode->Start = DivideInd+1;
	rNewNode->End = node->End;
	rNewNode->lptr = NULL;
	rNewNode->rptr = NULL;
	
	smallx = smally = smallz = HUGEDBL;
	bigx = bigy = bigz = -HUGEDBL;

	for(i=lNewNode->Start;i<=lNewNode->End;i++){
		if(rand[TreeInd[i]].x > bigx) bigx = rand[TreeInd[i]].x;
		if(rand[TreeInd[i]].x < smallx) smallx = rand[TreeInd[i]].x;
		if(rand[TreeInd[i]].y > bigy) bigy = rand[TreeInd[i]].y;
		if(rand[TreeInd[i]].y < smally) smally = rand[TreeInd[i]].y;
		if(rand[TreeInd[i]].z > bigz) bigz = rand[TreeInd[i]].z;
		if(rand[TreeInd[i]].z < smallz) smallz = rand[TreeInd[i]].z;
	}

	avgx = (smallx+bigx)/2.;
	avgy = (smally+bigy)/2.;
	avgz = (smallz+bigz)/2.;
	norm = sqrt(avgx*avgx + avgy*avgy + avgz*avgz);
	lNewNode->x = avgx/norm;
	lNewNode->y = avgy/norm;
	lNewNode->z = avgz/norm;

	mind = HUGEDBL;
	for(i=lNewNode->Start;i<=lNewNode->End;i++){
		d = lNewNode->x*rand[TreeInd[i]].x + lNewNode->y*rand[TreeInd[i]].y + lNewNode->z*rand[TreeInd[i]].z;
		if(d < mind) mind = d;
	}
	lNewNode->clowbound = mind;
	mindsq = mind*mind;
	lNewNode->c2lowbound = 2.*mindsq-1.;
	lNewNode->slowbound = sqrt(1-mindsq);

	x[0] = smallx;
	x[1] = bigx;
	y[0] = smally;
	y[1] = bigy;
	z[0] = smallz;
	z[1] = bigz;

	smallx = smally = smallz = HUGEDBL;
	bigx = bigy = bigz = -HUGEDBL;

	for(i=rNewNode->Start;i<=rNewNode->End;i++){
		if(rand[TreeInd[i]].x > bigx) bigx = rand[TreeInd[i]].x;
		if(rand[TreeInd[i]].x < smallx) smallx = rand[TreeInd[i]].x;
		if(rand[TreeInd[i]].y > bigy) bigy = rand[TreeInd[i]].y;
		if(rand[TreeInd[i]].y < smally) smally = rand[TreeInd[i]].y;
		if(rand[TreeInd[i]].z > bigz) bigz = rand[TreeInd[i]].z;
		if(rand[TreeInd[i]].z < smallz) smallz = rand[TreeInd[i]].z;
	}

	avgx = (smallx+bigx)/2.;
	avgy = (smally+bigy)/2.;
	avgz = (smallz+bigz)/2.;
	norm = sqrt(avgx*avgx + avgy*avgy + avgz*avgz);
	rNewNode->x = avgx/norm;
	rNewNode->y = avgy/norm;
	rNewNode->z = avgz/norm;

	mind = HUGEDBL;
	for(i=rNewNode->Start;i<=rNewNode->End;i++){
		d = rNewNode->x*rand[TreeInd[i]].x + rNewNode->y*rand[TreeInd[i]].y + rNewNode->z*rand[TreeInd[i]].z;
		if(d < mind) mind = d;
	}
	rNewNode->clowbound = mind;
	mindsq = mind*mind;
	rNewNode->c2lowbound = 2.*mindsq-1.;
	rNewNode->slowbound = sqrt(1-mindsq);

	angTreeBuildRan(lNewNode,x,y,z,rand,TreeInd,LogSamples,dsplit,Depth+1,dataname,out);
	x[0] = smallx;
	x[1] = bigx;
	y[0] = smally;
	y[1] = bigy;
	z[0] = smallz;
	z[1] = bigz;
	angTreeBuildRan(rNewNode,x,y,z,rand,TreeInd,LogSamples,dsplit,Depth+1,dataname,out);

	free(lNewNode);
	free(rNewNode);
	
	return;
}

int angTree(char filename[],splitlog dsplit[],int NumData,int NumSamples,int flag){

/* Coordinates tree building process, printing the results, and freeing memory allocated for a tree.	*
 * Also reorganizes data based on the tree and prints sample counts in the ".bin" data file.		*/

	int i,j,*TreeInd,LogSamples;
	char dfile[BUFFER_SIZE];
	double x[2],y[2],z[2];
	double avgx,avgy,avgz,norm,mind,mindsq,d;
	angTreeNode *root;
	pcsource *data;
	FILE *inout;
	#ifdef USE_DISK
	int pin;
	#endif

	if(NumSamples != 0){
		LogSamples = (int)(log((double)NumSamples)/M_LN2+0.5);
	}else{
		LogSamples = 0;
	}

	root = malloc(sizeof(angTreeNode));
	if(root == NULL){
		fprintf(stderr,"Failed to allocate root node\n");
		return 0;
	}

	#ifdef USE_DISK
	sprintf(dfile,"%s.bin.temp",filename);
	pin = open(dfile,O_RDONLY);
	data = mmap(NULL,NumData*sizeof(pcsource),PROT_READ,MAP_SHARED,pin,0);
	if((int)data == -1){
		fprintf(stderr,"Failed to map file to memory\n");
		return 0;
	}
	close(pin);
	#else
	sprintf(dfile,"%s0.bin",filename);
	inout = fopen(dfile,"r");
	fread(&i,sizeof(int),1,inout);
	fread(&NumData,sizeof(int),1,inout);
	for(i=0;i<NumSamples;i++){
		fread(&dsplit[i].Cnt,sizeof(int),1,inout);
	}

	data = malloc(NumData*sizeof(pcsource));
	if(data == NULL){
		fprintf(stderr,"Failled to allocate data in tree\n");
		return 0;
	}
	fread(data,sizeof(pcsource),NumData,inout);
	fclose(inout);
	remove(dfile);
	#endif
	
	TreeInd = malloc(NumData*sizeof(int));
	if(TreeInd == NULL){
		fprintf(stderr,"Failed to allocate TreeInd\n");
		return 0;
	}

	root->Sample = -1;
	root->Cnt = NumData;
	root->Start = 0;
	root->End = NumData-1;
	root->lptr = NULL;
	root->rptr = NULL;

	x[0] = y[0] = z[0] = HUGEDBL;
	x[1] = y[1] = z[1] = -HUGEDBL;

	for(i=0;i<NumData;i++){
		TreeInd[i] = i;
		if(data[i].x > x[1]) x[1] = data[i].x;
		if(data[i].x < x[0]) x[0] = data[i].x;
		if(data[i].y > y[1]) y[1] = data[i].y;
		if(data[i].y < y[0]) y[0] = data[i].y;
		if(data[i].z > z[1]) z[1] = data[i].z;
		if(data[i].z < z[0]) z[0] = data[i].z;
	}

	avgx = (x[0]+x[1])/2.;
	avgy = (y[0]+y[1])/2.;
	avgz = (z[0]+z[1])/2.;
	norm = sqrt(avgx*avgx + avgy*avgy + avgz*avgz);
	root->x = avgx/norm;
	root->y = avgy/norm;
	root->z = avgz/norm;

	mind = HUGEDBL;
	for(i=0;i<NumData;i++){
		d = root->x*data[i].x + root->y*data[i].y + root->z*data[i].z;
		if(d < mind) mind = d;
	}
	root->clowbound = mind;
	mindsq = mind*mind;
	root->c2lowbound = 2*mindsq-1.;
	root->slowbound = sqrt(1-mindsq);

	isplit = 0;
	CurSamp = 1;
	CurFile = 0;
	if(flag == 0){
		angTreeBuild(root,x,y,z,data,TreeInd,0,LogSamples,dsplit,filename,NULL);
	}else{
		angTreeBuildRan(root,x,y,z,data,TreeInd,LogSamples,dsplit,0,filename,NULL);
	}
	
	if(CurFile == 1){
		rewind(hlist->out);
		fwrite(&(hlist->Cnt),sizeof(int),1,hlist->out);
		fclose(hlist->out);
		if(root != NULL) free(root);
		sprintf(dfile,"%s0.bin",filename);
		inout = fopen(dfile,"wb");
		fwrite(&NumSamples,sizeof(int),1,inout);
		fwrite(&NumData,sizeof(int),1,inout);
		for(i=0;i<NumSamples;i++){
			fwrite(&dsplit[i].Cnt,sizeof(int),1,inout);
		}
		for(i=0;i<NumData;i++){
			fwrite(&(data[TreeInd[i]].x),sizeof(pcsource),1,inout);
		}
		fclose(inout);
		free(hlist);
	}else{
		if(root != NULL) free(root);
		for(i=0;i<CurFile;i++){
			rewind(hlist->out);
			fwrite(&(hlist->Cnt),sizeof(int),1,hlist->out);
			fclose(hlist->out);
			sprintf(dfile,"%s%d.bin",filename,i);
			inout = fopen(dfile,"wb");
			NumData = hlist->End - hlist->Start + 1;
			fwrite(&NumSamples,sizeof(int),1,inout);
			fwrite(&NumData,sizeof(int),1,inout);
			for(j=0;j<NumSamples;j++){
				fwrite(&dsplit[j].Cnt,sizeof(int),1,inout);
			}
			for(j=hlist->Start;j<=hlist->End;j++){
				fwrite(&(data[TreeInd[j]].x),sizeof(pcsource),1,inout);
			}
			fclose(inout);
			files = hlist->next;
			free(hlist);
			hlist = files;
		}
	}	

	free(TreeInd);
	
	#ifdef USE_DISK
	munmap(data,NumData*sizeof(pcsource));
	sprintf(dfile,"%s.bin.temp",filename);
	remove(dfile);
	#else
	fflush(stdout);
	free(data);
	#endif

	return CurFile;
}
