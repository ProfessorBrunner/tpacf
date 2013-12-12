#include "correlation.h"

/************************************************************************************

 Contains routines to allocate space for, read in, and free angular and spatial trees

 ************************************************************************************/


/***** Utility routines for both parallel and serial code *****/

static int getMaxTreesize(char *filelist){

/* Return size of largest tree in data sets listed in filelist */

        char dataname[BUFFER_SIZE],filename[BUFFER_SIZE];
        int i,numfiles,size,max;
        FILE *in,*list;

        #ifdef USE_MPI
        if(MyRank == 0){
        #endif
                
		max = 0;
                list = fopen(filelist,"r");
                if(list == NULL){
                        fprintf(stderr,"Failed to open list file in getMaxTreesize\n");
                        exit(1);
                }

                while(!feof(list)){
                        fscanf(list,"%s %d",&dataname[0],&numfiles);
                        if(!feof(list)){
                                for(i=0;i<numfiles;i++){
                                        sprintf(filename,"%s%d.tree",dataname,i);
                                        in = fopen(filename,"rb");
                                        if(in == NULL){
                                                fprintf(stderr,"Failed to open %s in getMaxTreesize\n",filename);
                                                exit(1);
                                        }
                                        fread(&size,sizeof(int),1,in);
                                        if(size>max) max = size;
                                        fclose(in);
                                }
                        }
                }

		fclose(list);
        #ifdef USE_MPI
        }
        MPI_Bcast(&max,1,MPI_INT,0,MPI_COMM_WORLD);
        #endif

	return max;
}

void treeGetMemory(void **block1, void **block2, char *filelist){

	int max,nodesize;
	
	max = getMaxTreesize(filelist);
	
	if(AngOrSpa == 0){
		nodesize = sizeof(angTreeNode);
	}else{
		nodesize = sizeof(spaTreeNode);
	}
	
	*block1 = malloc(max*nodesize);
	*block2 = malloc(max*nodesize);
	
	if(*block1 == NULL || *block2 == NULL){
		fprintf(stderr,"Failed to allocate tree space\n");
		exit(1);
	}
	
	return;
}


static int NodePos;	/* Used only in tree fill routines */

static spaTreeNode *spaTreeFill(spaTreeNode *root){

/* Fill in pointers in newly read tree */
	
	spaTreeNode *node;

	node = root+NodePos;

	if((int)node->lptr == 1){
		NodePos++;
		node->lptr = spaTreeFill(root);
		NodePos++;
		node->rptr = spaTreeFill(root);
	}else{
		node->lptr = NULL;
		node->rptr = NULL;
	}

	return node;
}

static angTreeNode *angTreeFill(angTreeNode *root){

/* Fill in pointers in newly read tree */

	angTreeNode *node;

	node = root+NodePos;

	if((int)node->lptr == 1){
		NodePos++;
		node->lptr = angTreeFill(root);
		NodePos++;
		node->rptr = angTreeFill(root);
	}else{
		node->lptr = NULL;
		node->rptr = NULL;
	}

	return node;
}

/************ Non-MPI only tree functions *******************/

#ifndef USE_MPI
void *treeRead(char dataName[], void *tblock){

/* Read in tree from file and fill pointers */

	int size;
	char treeFile[BUFFER_SIZE];
	FILE *in;
	
	TIMESTART(startTime);
	DEBUGPRINT("Entering treeRead\n");
	
	/* Append the dataName */
	sprintf(treeFile,"%s.tree",dataName);
	
	/* Open the file */
	if((in = fopen(treeFile,"r")) == NULL){
		fprintf(stderr,"Failed to open file %s\n",treeFile);
		exit(1);
	}

	/* Read in the size of the tree */
	if(fread(&size,sizeof(int),1,in) != 1){
		fprintf(stderr,"Failed to read header from %s\n",treeFile);
		exit(1);
	}
	
	/* Check if tblock points somewhere valid (sort of) */
	if(tblock == NULL){
		fprintf(stderr,"tree memory corrupted...aborting\n");
		exit(1);
	}
	
	/* Initialize NodePos for tree fill routines */
	NodePos = 0;
	
	/* Now read in the tree and fill in the pointers */
	if(AngOrSpa == 0){
		if(fread(tblock,sizeof(angTreeNode),size,in) != size){
			fprintf(stderr,"Failed to read tree from %s\n",treeFile);
			exit(1);
		}
		tblock = (void *)angTreeFill((angTreeNode *)tblock);
	}else{
		if(fread(tblock,sizeof(spaTreeNode),size,in) != size){
			fprintf(stderr,"Failed to read tree from %s\n",treeFile);
			exit(1);
		}
		tblock = (void *)spaTreeFill((spaTreeNode *)tblock);
	}
	
	fclose(in);
	
	DEBUGPRINT("Exiting treeRead\n");
	TIMESTOP(wallTreeRead,startTime);
	
	return tblock;
}
#endif


/* MPI-only tree functions */
 
#ifdef USE_MPI

static int spaTreeBuildMPI(MPI_Datatype *MPI_spaTreeNode){

/*  Builds MPI derived datatype for C datatype spaTreeNodedat, which is just	*
 *  spaTreeNode without the pointers.						*/

	int BlockLengths[12],ptrSize;
	MPI_Datatype typelist[12];
	MPI_Aint StartAddress,Address;
	MPI_Aint Displacements[12];
	spaTreeNode protonode;

	BlockLengths[0] = BlockLengths[1] = BlockLengths[2] = BlockLengths[3] = 1;
	BlockLengths[4] = BlockLengths[5] = BlockLengths[6] = BlockLengths[7] = 1;
	BlockLengths[8] = BlockLengths[9] = BlockLengths[10] = BlockLengths[11] = 1;
	typelist[0] = typelist[1] = typelist[2] = typelist[3] = MPI_INT;
	typelist[4] = typelist[5] = typelist[6] = MPI_DOUBLE;
	typelist[7] = typelist[8] = typelist[9] = MPI_DOUBLE;
	
	ptrSize = sizeof(spaTreeNode *);
	
	if(ptrSize == sizeof(int)){
		typelist[10] = typelist[11] = MPI_INT;
	}else if(ptrSize == sizeof(double)){
		typelist[10] = typelist[11] = MPI_DOUBLE;
	}else{
		printf("Cannot determine size of ptr to build MPI_spaTreeNode\n");
		return 0;
	}

	Displacements[0] = 0;
	MPI_Address(&(protonode.Sample),&StartAddress);
	MPI_Address(&(protonode.Cnt),&Address);
	Displacements[1] = Address - StartAddress;
	MPI_Address(&(protonode.Start),&Address);
	Displacements[2] = Address - StartAddress;
	MPI_Address(&(protonode.End),&Address);
	Displacements[3] = Address - StartAddress;
	MPI_Address(&(protonode.xmin),&Address);
	Displacements[4] = Address - StartAddress;
	MPI_Address(&(protonode.xmax),&Address);
	Displacements[5] = Address - StartAddress;
	MPI_Address(&(protonode.ymin),&Address);
	Displacements[6] = Address - StartAddress;
	MPI_Address(&(protonode.ymax),&Address);
	Displacements[7] = Address - StartAddress;
	MPI_Address(&(protonode.zmin),&Address);
	Displacements[8] = Address - StartAddress;
	MPI_Address(&(protonode.zmax),&Address);
	Displacements[9] = Address - StartAddress;
	MPI_Address(&(protonode.lptr),&Address);
	Displacements[10] = Address - StartAddress;
	MPI_Address(&(protonode.rptr),&Address);
	Displacements[11] = Address - StartAddress;
	
	MPI_Type_struct(12,BlockLengths,Displacements,typelist,MPI_spaTreeNode);
	MPI_Type_commit(MPI_spaTreeNode);
	
	return 1;
}

static int angTreeBuildMPI(MPI_Datatype *MPI_angTreeNode){

/*  Builds MPI derived datatype for C datatype angTreeNode */

	int BlockLengths[12],ptrSize;
	MPI_Datatype typelist[12];
	MPI_Aint StartAddress,Address;
	MPI_Aint Displacements[12];
	angTreeNode protonode;

	BlockLengths[0] = BlockLengths[1] = BlockLengths[2] = BlockLengths[3] = 1;
	BlockLengths[4] = BlockLengths[5] = BlockLengths[6] = BlockLengths[7] = 1;
	BlockLengths[8] = BlockLengths[9] = BlockLengths[10] = BlockLengths[11] = 1;
	typelist[0] = typelist[1] = typelist[2] = typelist[3] = MPI_INT;
	typelist[4] = typelist[5] = typelist[6] = MPI_DOUBLE;
	typelist[7] = typelist[8] = typelist[9] = MPI_DOUBLE;

	ptrSize = sizeof(angTreeNode *);
	
	
	if(ptrSize == sizeof(int)){
		typelist[10] = typelist[11] = MPI_INT;
	}else if(ptrSize == sizeof(double)){
		typelist[10] = typelist[11] = MPI_DOUBLE;
	}else{
		printf("Cannot determine size of ptr to build MPI_angTreeNode\n");
		return 0;
	}

	Displacements[0] = 0;
	MPI_Address(&(protonode.Sample),&StartAddress);
	MPI_Address(&(protonode.Cnt),&Address);
	Displacements[1] = Address - StartAddress;
	MPI_Address(&(protonode.Start),&Address);
	Displacements[2] = Address - StartAddress;
	MPI_Address(&(protonode.End),&Address);
	Displacements[3] = Address - StartAddress;
	MPI_Address(&(protonode.x),&Address);
	Displacements[4] = Address - StartAddress;
	MPI_Address(&(protonode.y),&Address);
	Displacements[5] = Address - StartAddress;
	MPI_Address(&(protonode.z),&Address);
	Displacements[6] = Address - StartAddress;
	MPI_Address(&(protonode.clowbound),&Address);
	Displacements[7] = Address - StartAddress;
	MPI_Address(&(protonode.c2lowbound),&Address);
	Displacements[8] = Address - StartAddress;
	MPI_Address(&(protonode.slowbound),&Address);
	Displacements[9] = Address - StartAddress;
	MPI_Address(&(protonode.lptr),&Address);
	Displacements[10] = Address - StartAddress;
	MPI_Address(&(protonode.rptr),&Address);
	Displacements[11] = Address - StartAddress;
	
	MPI_Type_struct(12,BlockLengths,Displacements,typelist,MPI_angTreeNode);
	MPI_Type_commit(MPI_angTreeNode);
	
	return 1;
}

void *treeRead(char dataName[], void *tblock){

/* Reads in tree to all MPI processes */

	static int firstc = 1;
	static MPI_Datatype MPI_spaTreeNode;
	static MPI_Datatype MPI_angTreeNode;
	int size, ierr;
	char treeFile[BUFFER_SIZE];
	MPI_File in;
	MPI_Status status;

	TIMESTART(startTime);
	
	DEBUGPRINT("Entering treeRead\n");

	/* If this is the first time through, build the MPI datatype for a tree node */
	if(firstc){
		if(AngOrSpa == 0){
			if(!angTreeBuildMPI(&MPI_angTreeNode)){
				fprintf(stderr,"Failed to build MPI_angTreeNode\n");
				exit(1);
			}
		}else{
			if(!spaTreeBuildMPI(&MPI_spaTreeNode)){
				fprintf(stderr,"Failed to build MPI_spaTreeNode\n");
				exit(1);
			}
		}
		firstc = 0;
	}
	
	/* Append the dataName */
	sprintf(treeFile,"%s.tree",dataName);
	
	/* Open the file */
	if((ierr = MPI_File_open(MPI_COMM_WORLD,treeFile,MPI_MODE_RDONLY,MPI_INFO_NULL,&in)) != MPI_SUCCESS){
		fprintf(stderr,"Failed to open tree file %s\n",treeFile);
		exit(1);
	}

	/* Read in the tree size */
	if((ierr = MPI_File_read_all(in,&size,1,MPI_INT,&status)) != MPI_SUCCESS){
		fprintf(stderr,"Failed to read tree size from %s\n",treeFile);
		exit(1);
	}

	/* Read in the tree */
	if(AngOrSpa == 0){
		ierr = MPI_File_read_all(in,tblock,size,MPI_angTreeNode,&status);
	}else{
		ierr = MPI_File_read_all(in,tblock,size,MPI_spaTreeNode,&status);
	}
	
	/* Make sure the read succeeded */
	if(ierr != MPI_SUCCESS){
		fprintf(stderr,"Failed to read tree from %s\n",treeFile);
		exit(1);
	}
	
	MPI_File_close(&in);	

	/* Initialize NodePos for tree fill routines */
	NodePos = 0;
	
	/* Now fill in the pointers */
	if(AngOrSpa == 0){
		tblock = (void *)angTreeFill((angTreeNode *)tblock);
	}else{
		tblock = (void *)spaTreeFill((spaTreeNode *)tblock);
	}

	DEBUGPRINT("Exiting treeRead\n");
	
	TIMESTOP(wallTreeRead,startTime);

	return tblock;
}
#endif


/************ Parallel (MPI or OpenMP) tree functions ************/

#if defined (USE_MPI) || defined (USE_OMP)
int spaTreeWorkNodes(spaTreeNode *node,int WorkLevel,int CurNum,spaTreeNode *worknodes[]){

/*  Gets work nodes in tree for parallel distribution.  Fills worknodes	*
 *  with the nodes at depth WORKLEVEL					*/

	int Depth;
	typedef struct stack_{
		spaTreeNode *node;
		int Depth;
		struct stack_ *next;
	}stack;
	stack *mystack,*newstack;
	
	CurNum = 0;
	
	mystack = malloc(sizeof(stack));
	if(mystack == NULL){
		printf("Failed to allocate stack in spaTreeWorkNodes...\n");
		return 0;
	}
	mystack->node = node;
	mystack->Depth = 0;
	mystack->next = NULL;
	
	while(mystack != NULL){
	
		node = mystack->node;
		Depth = mystack->Depth;

		if(node == NULL){
			newstack = mystack->next;
			free(mystack);
			mystack = newstack;
			continue;
		}

		if(Depth < WorkLevel){
			newstack = malloc(sizeof(stack));
			if(newstack == NULL){
				printf("Failed to allocate stack in spaTreeWorkNodes\n");
				return 0;
			}
			newstack->node = node->rptr;
			newstack->Depth = Depth + 1;
			newstack->next = mystack->next;
			mystack->node = node->lptr;
			mystack->Depth = Depth + 1;
			mystack->next = newstack;
		}else if(Depth == WorkLevel){
			worknodes[CurNum] = node;
			CurNum++;
			newstack = mystack->next;
			free(mystack);
			mystack = newstack;
		}
	}
	return 1;
}

int angTreeWorkNodes(angTreeNode *node,int WorkLevel,int CurNum,angTreeNode *worknodes[]){

/*  Gets work nodes in tree for parallel distribution.  Fills worknodes	*
 *  with the nodes at depth WORKLEVEL					*/

	int Depth;
	typedef struct stack_{
		angTreeNode *node;
		int Depth;
		struct stack_ *next;
	}stack;
	stack *mystack,*newstack;
	
	CurNum = 0;
	
	mystack = malloc(sizeof(stack));
	if(mystack == NULL){
		printf("Failed to allocate stack in angTreeWorkNodes...\n");
		return 0;
	}
	mystack->node = node;
	mystack->Depth = 0;
	mystack->next = NULL;
	
	while(mystack != NULL){
	
		node = mystack->node;
		Depth = mystack->Depth;

		if(node == NULL){
			newstack = mystack->next;
			free(mystack);
			mystack = newstack;
			continue;
		}

		if(Depth < WorkLevel){
			newstack = malloc(sizeof(stack));
			if(newstack == NULL){
				printf("Failed to allocate stack in angTreeWorkNodes\n");
				return 0;
			}
			newstack->node = node->rptr;
			newstack->Depth = Depth + 1;
			newstack->next = mystack->next;
			mystack->node = node->lptr;
			mystack->Depth = Depth + 1;
			mystack->next = newstack;
		}else if(Depth == WorkLevel){
			worknodes[CurNum] = node;
			CurNum++;
			newstack = mystack->next;
			free(mystack);
			mystack = newstack;
		}
	}
	return 1;
}
#endif
