#include "correlation.h"

#define FALSE 0


/* Some routines to get bin counts using the sure-fire brute force algorithm */
/*
static void bruteACspatial(const pcsource data[],bin bins[],int NumData){

	int i,j,Bin1;
	double dx,dy,dz,dist;

	for(i=0;i<NumData;i++){
		for(j=i+1;j<NumData;j++){
			dx = data[i].x - data[j].x;
			dy = data[i].y - data[j].y;
			dz = data[i].z - data[j].z;
			dist = dx*dx + dy*dy + dz*dz;
			Bin1 = 0;
			while(dist < bins[Bin1].limit){
				Bin1++;
			}
			bins[Bin1].Cnt[0]++;
		}
	}
	
	return;
}

static void bruteCCspatial(const pcsource data[],const pcsource rand[],bin bins[],int NumData,int NumRand){

	int i,j,Bin1;
	double dx,dy,dz,dist;
	
	for(i=0;i<NumData;i++){
		for(j=0;j<NumRand;j++){
			dx = data[i].x - rand[j].x;
			dy = data[i].y - rand[j].y;
			dz = data[i].z - rand[j].z;
			dist = dx*dx + dy*dy + dz*dz;
			Bin1 = 0;
			while(dist < bins[Bin1].limit){
				Bin1++;
			}
			bins[Bin1].Cnt[0]++;
		}
	}
	
	return;
}

static void bruteACangular(const pcsource data[],bin bins[],int NumData){

	int i,j,Bin1;
	double csep;

	for(i=0;i<NumData;i++){
		for(j=i+1;j<NumData;j++){
			csep = data[i].x*data[j].x + data[i].y*data[j].y + data[i].z*data[j].z;
			Bin1 = 0;
			while(csep > bins[Bin1].limit){
				Bin1++;
			}
			bins[Bin1].Cnt[0]++;
		}
	}
	
	return;
}

static void bruteCCangular(const pcsource data[],const pcsource rand[],bin bins[],int NumData,int NumRand){

	int i,j,Bin1;
	double csep;

	for(i=0;i<NumData;i++){
		for(j=0;j<NumRand;j++){
			csep = data[i].x*rand[j].x + data[i].y*rand[j].y + data[i].z*rand[j].z;
			Bin1 = 0;
			while(csep > bins[Bin1].limit){
				Bin1++;
			}
			bins[Bin1].Cnt[0]++;
		}
	}

	return;
}
*/
/* End SLOW brute force routines */


/****************************************************************************************
 *											*
 *				dualtreeACspatial					*
 *											*
 ****************************************************************************************/
static void dualtreeACspatial(const pcsource data[],void *nodeA,void *nodeB,bin bins[]){

/* Main routine to autocorrelate a data set.	*
 * Described in Dolence & Brunner (2007).	*/

	unsigned int i,j,N2Start,N2End,Sample1,Sample2,BinEnd;
	unsigned long long int Add,*CntSave;
	double mindist,maxdist,dxmin,dymin,dzmin,dxmax,dymax,dzmax;
	int Bin1,Bin2,BinStart,BinStop;
	spaTreeNode *node1,*node2;
	typedef struct stack_{
		spaTreeNode *n1,*n2;
		int Fbin,Lbin;
		struct stack_ *next;
	}stack;
	stack *mystack,*newstack;

	BinStart = 0;
	BinStop = NumBins;

	node1 = (spaTreeNode *)nodeA;
	node2 = (spaTreeNode *)nodeB;
	
	CntSave = malloc((BinStop+2)*sizeof(unsigned long long int));
	if(CntSave == NULL){
		printf("Failed to allocate CntSave in dualtreeAC...WRONG RESULTS\n");
		return;
	}
	
	mystack = malloc(sizeof(stack));
	if(mystack == NULL){
		printf("Failed to allocate stack in dualtreeAC...WRONG RESULTS\n");
		return;
	}
	mystack->n1 = node1;
	mystack->n2 = node2;
	mystack->Fbin = BinStart;
	mystack->Lbin = BinStop;
	mystack->next = NULL;
	
	while(mystack != NULL){	/* Loop until all pairs have been accounted for */
	
		node1 = mystack->n1;
		node2 = mystack->n2;
		BinStart = mystack->Fbin;
		BinStop = mystack->Lbin;

		/* Avoid redundant calculations */
		if(node1->Start >= node2->End){
			newstack = mystack->next;
			free(mystack);
			mystack = newstack;
			continue;
		}

		/* Special case if node1 and node2 are the same */
		if(node1 == node2){
			if(node1->lptr == NULL){	/* node is a leaf, must do slow brute force calculation */
				/* Save the sample counts */
				for(i=BinStart;i<=BinStop;i++){
					CntSave[i] = bins[i].Cnt[0];
				}
				
				/* More efficient to avoid references to a node in the inner loop */
				N2End = node1->End;
				for(i=node1->Start;i<=node1->End;i++){
					for(j=i+1;j<=N2End;j++){
						dxmin = data[i].x - data[j].x;
						dymin = data[i].y - data[j].y;
						dzmin = data[i].z - data[j].z;
						mindist = dxmin*dxmin + dymin*dymin + dzmin*dzmin;
						Bin1 = BinStart;
						while(mindist < bins[Bin1].limit){
							Bin1++;
						}
						bins[Bin1].Cnt[0]++;
					}
				}
				Sample1 = node1->Sample;
				for(i=BinStart;i<=BinStop;i++){	/* Update sample counts */
					bins[i].Cnt[Sample1] += (bins[i].Cnt[0] - CntSave[i]);
				}
				/* Move on to the next set of nodes in the stack */
				newstack = mystack->next;
				free(mystack);
				mystack = newstack;
			}else{				/* node is not a leaf, so explore its children */
				newstack = malloc(sizeof(stack));
				if(newstack == NULL){
					printf("Failed to allocate stack in dualtreeAC...WRONG RESULTS\n");
					return;
				}
				newstack->n1 = node1->rptr;
				newstack->n2 = node2;
				newstack->Fbin = BinStart;
				newstack->Lbin = BinStop;
				newstack->next = mystack->next;
				mystack->n1 = node1->lptr;
				mystack->n2 = node2;
				mystack->Fbin = BinStart;
				mystack->Lbin = BinStop;
				mystack->next = newstack;
			}
			continue;
		}

		if(node1->xmax < node2->xmin){
			dxmin = node2->xmin - node1->xmax;;
			dxmax = node2->xmax - node1->xmin;
		}else if(node1->xmin > node2->xmax){
			dxmin = node1->xmin - node2->xmax;
			dxmax = node1->xmax - node2->xmin;
		}else{
			dxmin = node1->xmax - node2->xmin;
			dxmax = node2->xmax - node1->xmin;
			if(dxmin > dxmax) dxmax = dxmin;
			dxmin = 0.;
		}

		if(node1->ymax < node2->ymin){
			dymin = node2->ymin - node1->ymax;
			dymax = node2->ymax - node1->ymin;
		}else if(node1->ymin > node2->ymax){
			dymin = node1->ymin - node2->ymax;
			dymax = node1->ymax - node2->ymin;
		}else{
			dymin = node1->ymax - node2->ymin;
			dymax = node2->ymax - node1->ymin;
			if(dymin > dymax) dymax = dymin;
			dymin = 0.;
		}	

		if(node1->zmax < node2->zmin){
			dzmin = node2->zmin - node1->zmax;
			dzmax = node2->zmax - node1->zmin;
		}else if(node1->zmin > node2->zmax){
			dzmin = node1->zmin - node2->zmax;
			dzmax = node1->zmax - node2->zmin;
		}else{
			dzmin = node1->zmax - node2->zmin;
			dzmax = node2->zmax - node1->zmin;
			if(dzmin > dzmax) dzmax = dzmin;
			dzmin = 0.;
		}

		mindist = dxmin*dxmin + dymin*dymin + dzmin*dzmin;
		
		/* check if nodes are too distant to be of interest */
		if(mindist > bins[0].limit){
			newstack = mystack->next;
			free(mystack);
			mystack = newstack;
			continue;
		}

		/* Bin1 is distance bin for smallest separation */
		Bin1 = BinStart;
		while(mindist < bins[Bin1].limit){
			Bin1++;
		}

		maxdist = dxmax*dxmax + dymax*dymax + dzmax*dzmax;

		/* Bin2 is distance bin for largest separation */
		Bin2 = BinStart;
		while(maxdist < bins[Bin2].limit){
			Bin2++;
		}


		if(Bin1 == Bin2){	/* Bins are the same, so add node1.Cnt*node2.Cnt to total counts and sample counts */
			Add = node1->Cnt*node2->Cnt;
			bins[Bin1].Cnt[0] += Add;
			Sample1 = node1->Sample;
			Sample2 = node2->Sample;
			if(Sample1 == Sample2){
				bins[Bin1].Cnt[Sample1] += Add;
			}else{
				bins[Bin1].Cnt[Sample1] += Add;
				bins[Bin1].Cnt[Sample2] += Add;
			}
		
			/* Move on to the next set of nodes in the stack */
			newstack = mystack->next;
			free(mystack);
			mystack = newstack;
		
		}else if(node1->lptr == NULL && node2->lptr == NULL){	/* nodes are leafs, try single tree subsumption, then brute force */

			BinStart = Bin2;	/* Most distant bin possible */
			BinEnd = Bin1;		/* Closest bin possible */
			
			for(i=BinStart;i<=BinEnd;i++){	/* Save sample counts */
				CntSave[i] = bins[i].Cnt[0];
			}
			
			N2Start = node2->Start;
			N2End = node2->End;
			for(i=node1->Start;i<=node1->End;i++){
			
				if(data[i].x < node2->xmin){
					dxmin = node2->xmin - data[i].x;
					dxmax = node2->xmax - data[i].x;
				}else if(data[i].x > node2->xmax){
					dxmin = data[i].x - node2->xmax;
					dxmax = data[i].x - node2->xmin;
				}else{
					dxmin = data[i].x - node2->xmin;
					dxmax = node2->xmax - data[i].x;
					if(dxmin > dxmax) dxmax = dxmin;
					dxmin = 0.;
				}
				if(data[i].y < node2->ymin){
					dymin = node2->ymin - data[i].y;
					dymax = node2->ymax - data[i].y;
				}else if(data[i].y > node2->ymax){
					dymin = data[i].y - node2->ymax;
					dymax = data[i].y - node2->ymin;
				}else{
					dymin = data[i].y - node2->ymin;
					dymax = node2->ymax - data[i].y;
					if(dymin > dymax) dymax = dymin;
					dymin = 0.;
				}
				if(data[i].z < node2->zmin){
					dzmin = node2->zmin - data[i].z;
					dzmax = node2->zmax - data[i].z;
				}else if(data[i].z > node2->zmax){
					dzmin = data[i].z - node2->zmax;
					dzmax = data[i].z - node2->zmin;
				}else{
					dzmin = data[i].z - node2->zmin;
					dzmax = node2->zmax - data[i].z;
					if(dzmin > dzmax) dzmax = dzmin;
					dzmin = 0.;
				}

				maxdist = dxmax*dxmax + dymax*dymax + dzmax*dzmax;
				mindist = dxmin*dxmin + dymin*dymin + dzmin*dzmin;
			
				Bin1 = BinStart;
				while(mindist < bins[Bin1].limit){
					Bin1++;
				}
				Bin2 = BinStart;
				while(maxdist < bins[Bin2].limit){
					Bin2++;
				}

				/* Bins are the same, so avoid looking at the individual points in node2 by	*
				 * single-tree subsumption							*/
				if(Bin1 == Bin2){
					bins[Bin1].Cnt[0] += node2->Cnt;
					continue;
				}
				
				/* Brute force calculation necessary */
				/* More efficient to avoid references to a node in the inner loop */
				for(j=N2Start;j<=N2End;j++){
					dxmin = data[i].x - data[j].x;
					dymin = data[i].y - data[j].y;
					dzmin = data[i].z - data[j].z;
					mindist = dxmin*dxmin + dymin*dymin + dzmin*dzmin;
					Bin1 = Bin2;
					while(mindist < bins[Bin1].limit){
						Bin1++;
					}
					bins[Bin1].Cnt[0]++;
				}
			}
			
			/* Update sample counts */
			Sample1 = node1->Sample;
			Sample2 = node2->Sample;
			if(Sample1 == Sample2){
				for(i=BinStart;i<=BinEnd;i++){
					bins[i].Cnt[Sample1] += bins[i].Cnt[0] - CntSave[i];
				}
			}else{
				for(i=BinStart;i<=BinEnd;i++){
					Add = bins[i].Cnt[0] - CntSave[i];
					bins[i].Cnt[Sample1] += Add;
					bins[i].Cnt[Sample2] += Add;
				}
			}
			
			/* Move on to the next set of nodes in the stack */
			newstack = mystack->next;
			free(mystack);
			mystack = newstack;
		
		}else{	/* Explore children */
		
			if(node1->Cnt > node2->Cnt){	/* open node1 */
				newstack = malloc(sizeof(stack));
				if(newstack == NULL){
					printf("Failed to allocate stack in dualtreeAC...WRONG RESULTS\n");
					return;
				}
				newstack->n1 = node1->rptr;
				newstack->n2 = node2;
				newstack->Fbin = Bin2;
				newstack->Lbin = Bin1;
				newstack->next = mystack->next;
				mystack->n1 = node1->lptr;
				mystack->n2 = node2;
				mystack->Fbin = Bin2;
				mystack->Lbin = Bin1;
				mystack->next = newstack;

			}else{	/* open node2 */
				newstack = malloc(sizeof(stack));
				if(newstack == NULL){
					printf("Failed to allocate stack in dualtreeAC...WRONG RESULTS\n");
					return;
				}
				newstack->n1 = node1;
				newstack->n2 = node2->rptr;
				newstack->Fbin = Bin2;
				newstack->Lbin = Bin1;
				newstack->next = mystack->next;
				mystack->n1 = node1;
				mystack->n2 = node2->lptr;
				mystack->Fbin = Bin2;
				mystack->Lbin = Bin1;
				mystack->next = newstack;
			}
		}
	}

	free(CntSave);

	/* All done! */
	
	return;

}

/****************************************************************************************
 *											*
 *				      dualtreeCCspatial					*
 *											*
 ****************************************************************************************/
static void dualtreeCCspatial(const pcsource data1[],const pcsource data2[],void *nodeA,void *nodeB,bin bins[]){

/* Main routine to cross correlate two data sets A and B.	*
 * Described in Dolence & Brunner (2007).			*/

	unsigned int i,j,N2Start,N2End,Sample1,Sample2,BinEnd;
	unsigned long long int Add,*CntSave;
	double mindist,maxdist,dxmin,dymin,dzmin,dxmax,dymax,dzmax;
	int Bin1,Bin2,BinStart,BinStop;
	spaTreeNode *node1,*node2;
	typedef struct stack_{
		spaTreeNode *n1,*n2;
		int Fbin,Lbin;
		struct stack_ *next;
	}stack;
	stack *mystack,*newstack;
	
	BinStart = 0;
	BinStop = NumBins;

	node1 = (spaTreeNode *)nodeA;
	node2 = (spaTreeNode *)nodeB;

	CntSave = malloc((BinStop+2)*sizeof(unsigned long long int));
	if(CntSave == NULL){
		printf("Failed to allocate CntSave in dualtreeCC...WRONG RESULTS\n");
		return;
	}
	
	mystack = malloc(sizeof(stack));
	if(mystack == NULL){
		printf("Failed to allocate stack in dualtreeCC...WRONG RESULTS\n");
		return;
	}
	mystack->n1 = node1;
	mystack->n2 = node2;
	mystack->Fbin = BinStart;
	mystack->Lbin = BinStop;
	mystack->next = NULL;
	
	while(mystack != NULL){	/* Loop until all pairs have been accounted for */
	
		node1 = mystack->n1;
		node2 = mystack->n2;
		BinStart = mystack->Fbin;
		BinStop = mystack->Lbin;

		if(node1->xmax < node2->xmin){
			dxmin = node2->xmin - node1->xmax;;
			dxmax = node2->xmax - node1->xmin;
		}else if(node1->xmin > node2->xmax){
			dxmin = node1->xmin - node2->xmax;
			dxmax = node1->xmax - node2->xmin;
		}else{
			dxmin = node1->xmax - node2->xmin;
			dxmax = node2->xmax - node1->xmin;
			if(dxmin > dxmax)	dxmax = dxmin;
			dxmin = 0.;
		}

		if(node1->ymax < node2->ymin){
			dymin = node2->ymin - node1->ymax;
			dymax = node2->ymax - node1->ymin;
		}else if(node1->ymin > node2->ymax){
			dymin = node1->ymin - node2->ymax;
			dymax = node1->ymax - node2->ymin;
		}else{
			dymin = node1->ymax - node2->ymin;
			dymax = node2->ymax - node1->ymin;
			if(dymin > dymax) dymax = dymin;
			dymin = 0.;
		}	

		if(node1->zmax < node2->zmin){
			dzmin = node2->zmin - node1->zmax;
			dzmax = node2->zmax - node1->zmin;
		}else if(node1->zmin > node2->zmax){
			dzmin = node1->zmin - node2->zmax;
			dzmax = node1->zmax - node2->zmin;
		}else{
			dzmin = node1->zmax - node2->zmin;
			dzmax = node2->zmax - node1->zmin;
			if(dzmin > dzmax) dzmax = dzmin;
			dzmin = 0.;
		}

		mindist = dxmin*dxmin + dymin*dymin + dzmin*dzmin;
		
		/* check if nodes are too distant to be of interest */
		if(mindist > bins[0].limit){
			newstack = mystack->next;
			free(mystack);
			mystack = newstack;
			continue;
		}

		/* Bin1 is distance bin for smallest separation */
		Bin1 = BinStart;
		while(mindist < bins[Bin1].limit){
			Bin1++;
		}

		maxdist = dxmax*dxmax + dymax*dymax + dzmax*dzmax;

		/* Bin2 is distance bin for largest separation */
		Bin2 = BinStart;
		while(maxdist < bins[Bin2].limit){
			Bin2++;
		}


		if(Bin1 == Bin2){	/* Bins are the same, so add node1.Cnt*node2.Cnt to total counts and sample counts */
			Add = node1->Cnt*node2->Cnt;
			bins[Bin1].Cnt[0] += Add;
			Sample1 = node1->Sample;
			Sample2 = node2->Sample;
			if(Sample1 == Sample2){
				bins[Bin1].Cnt[Sample1] += Add;
			}else{
				bins[Bin1].Cnt[Sample1] += Add;
				bins[Bin1].Cnt[Sample2] += Add;
			}
		
			/* Move on to the next set of nodes in the stack */
			newstack = mystack->next;
			free(mystack);
			mystack = newstack;
		
		}else if(node1->lptr == NULL && node2->lptr == NULL){	/* nodes are leafs, try single tree subsumption, then brute force */

			BinStart = Bin2;	/* Most distant bin possible */
			BinEnd = Bin1;		/* Closest bin possible */
			
			for(i=BinStart;i<=BinEnd;i++){	/* Save sample counts */
				CntSave[i] = bins[i].Cnt[0];
			}

			N2Start = node2->Start;
			N2End = node2->End;
			for(i=node1->Start;i<=node1->End;i++){
			
				if(data1[i].x < node2->xmin){
					dxmin = node2->xmin - data1[i].x;
					dxmax = node2->xmax - data1[i].x;
				}else if(data1[i].x > node2->xmax){
					dxmin = data1[i].x - node2->xmax;
					dxmax = data1[i].x - node2->xmin;
				}else{
					dxmin = data1[i].x - node2->xmin;
					dxmax = node2->xmax - data1[i].x;
					if(dxmin > dxmax) dxmax = dxmin;
					dxmin = 0.;
				}
				if(data1[i].y < node2->ymin){
					dymin = node2->ymin - data1[i].y;
					dymax = node2->ymax - data1[i].y;
				}else if(data1[i].y > node2->ymax){
					dymin = data1[i].y - node2->ymax;
					dymax = data1[i].y - node2->ymin;
				}else{
					dymin = data1[i].y - node2->ymin;
					dymax = node2->ymax - data1[i].y;
				if(dymin > dymax) dymax = dymin;
					dymin = 0.;
				}
				if(data1[i].z < node2->zmin){
					dzmin = node2->zmin - data1[i].z;
					dzmax = node2->zmax - data1[i].z;
				}else if(data1[i].z > node2->zmax){
					dzmin = data1[i].z - node2->zmax;
					dzmax = data1[i].z - node2->zmin;
				}else{
					dzmin = data1[i].z - node2->zmin;
					dzmax = node2->zmax - data1[i].z;
					if(dzmin > dzmax) dzmax = dzmin;
					dzmin = 0.;
				}

				maxdist = dxmax*dxmax + dymax*dymax + dzmax*dzmax;
				mindist = dxmin*dxmin + dymin*dymin + dzmin*dzmin;
				
				Bin1 = BinStart;
				while(mindist < bins[Bin1].limit){
					Bin1++;
				}
				Bin2 = BinStart;
				while(maxdist < bins[Bin2].limit){
					Bin2++;
				}

				/* Bins are the same, so avoid looking at the individual points in node2 by	*
				 * single-tree subsumption							*/
				if(Bin1 == Bin2){
					bins[Bin1].Cnt[0] += node2->Cnt;
					continue;
				}
				
				/* Brute force calculation necessary */
				/* More efficient to avoid references to a node in the inner loop */
				for(j=N2Start;j<=N2End;j++){
					dxmin = data1[i].x - data2[j].x;
					dymin = data1[i].y - data2[j].y;
					dzmin = data1[i].z - data2[j].z;
					mindist = dxmin*dxmin + dymin*dymin + dzmin*dzmin;
					Bin1 = Bin2;
					while(mindist < bins[Bin1].limit){
						Bin1++;
					}
					bins[Bin1].Cnt[0]++;
				}
			}
			
			/* Update sample counts */
			Sample1 = node1->Sample;
			Sample2 = node2->Sample;
			if(Sample1 == Sample2){
				for(i=BinStart;i<=BinEnd;i++){
					bins[i].Cnt[Sample1] += bins[i].Cnt[0] - CntSave[i];
				}
			}else{
				for(i=BinStart;i<=BinEnd;i++){
					Add = bins[i].Cnt[0] - CntSave[i];
					bins[i].Cnt[Sample1] += Add;
					bins[i].Cnt[Sample2] += Add;
				}
			}
			
			/* Move on to the next set of nodes in the stack */
			newstack = mystack->next;
			free(mystack);
			mystack = newstack;
		
		}else{	/* Explore children */
		
			if(node1->Cnt > node2->Cnt){	/* open node1 */
				newstack = malloc(sizeof(stack));
				if(newstack == NULL){
					printf("Failed to allocate stack in dualtreeCC...WRONG RESULTS\n");
					return;
				}
				newstack->n1 = node1->rptr;
				newstack->n2 = node2;
				newstack->Fbin = Bin2;
				newstack->Lbin = Bin1;
				newstack->next = mystack->next;
				mystack->n1 = node1->lptr;
				mystack->n2 = node2;
				mystack->Fbin = Bin2;
				mystack->Lbin = Bin1;
				mystack->next = newstack;

			}else{	/* open node2 */
				newstack = malloc(sizeof(stack));
				if(newstack == NULL){
					printf("Failed to allocate stack in dualtreeCC...WRONG RESULTS\n");
					return;
				}
				newstack->n1 = node1;
				newstack->n2 = node2->rptr;
				newstack->Fbin = Bin2;
				newstack->Lbin = Bin1;
				newstack->next = mystack->next;
				mystack->n1 = node1;
				mystack->n2 = node2->lptr;
				mystack->Fbin = Bin2;
				mystack->Lbin = Bin1;
				mystack->next = newstack;
			}
		}
	}

	free(CntSave);

	/* All done! */
	
	return;

}

/* Dot product */
#define DOTP(a,b) (((a).x)*((b).x)+((a).y)*((b).y)+((a).z)*((b).z))

/****************************************************************************************
 *											*
 *					dualtreeACangular				*
 *											*
 ****************************************************************************************/
static void dualtreeACangular(const pcsource data[],void *nodeA,void *nodeB,bin bins[]){

/* Main routine to autocorrelate a data set.	*
 * Described in Dolence & Brunner (2007).	*/

	unsigned int i,j,N2Start,N2End,Sample1,Sample2,BinEnd;
	unsigned long long int Add,*CntSave;
	double mindist,maxdist,csep,ssep,crange,srange,term1,term2;
	int Bin1,Bin2,BinStart,BinStop;
	angTreeNode *node1,*node2;
	typedef struct stack_{
		angTreeNode *n1,*n2;
		int Fbin,Lbin;
		struct stack_ *next;
	}stack;
	stack *mystack,*newstack;

	BinStart = 0;
	BinStop = NumBins;

	node1 = (angTreeNode *)nodeA;
	node2 = (angTreeNode *)nodeB;
	
	CntSave = malloc((BinStop+2)*sizeof(unsigned long long int));
	if(CntSave == NULL){
		printf("Failed to allocate CntSave in dualtreeAC...WRONG RESULTS\n");
		return;
	}
	
	mystack = malloc(sizeof(stack));
	if(mystack == NULL){
		printf("Failed to allocate stack in dualtreeAC...WRONG RESULTS\n");
		return;
	}
	mystack->n1 = node1;
	mystack->n2 = node2;
	mystack->Fbin = BinStart;
	mystack->Lbin = BinStop;
	mystack->next = NULL;
	
	while(mystack != NULL){	/* Loop until all pairs have been accounted for */
	
		node1 = mystack->n1;
		node2 = mystack->n2;
		BinStart = mystack->Fbin;
		BinStop = mystack->Lbin;

		/* Avoid redundant calculations */
		if(node1->Start >= node2->End){
			newstack = mystack->next;
			free(mystack);
			mystack = newstack;
			continue;
		}

		/* Special case if node1 and node2 are the same */
		if(node1 == node2){
			Bin2 = BinStart;
			while(node1->c2lowbound > bins[Bin2].limit){	/* Figure out where we need to start the bin search */
				Bin2++;
			}
			if(node1->lptr == NULL){	/* node is a leaf, must do slow brute force calculation */
				/* Save the sample counts */
				for(i=Bin2;i<=BinStop;i++){
					CntSave[i] = bins[i].Cnt[0];
				}
				
				/* More efficient to avoid references to a node in the inner loop */
				N2End = node1->End;
				for(i=node1->Start;i<=node1->End;i++){
					for(j=i+1;j<=N2End;j++){
						csep = DOTP(data[i],data[j]);
						Bin1 = Bin2;
						while(csep > bins[Bin1].limit){
							Bin1++;
						}
						bins[Bin1].Cnt[0]++;
					}
				}
				Sample1 = node1->Sample;
				for(i=Bin2;i<=BinStop;i++){	/* Update sample counts */
					bins[i].Cnt[Sample1] += (bins[i].Cnt[0] - CntSave[i]);
				}
				/* Move on to the next set of nodes in the stack */
				newstack = mystack->next;
				free(mystack);
				mystack = newstack;
			}else{				/* node is not a leaf, so explore its children */
				newstack = malloc(sizeof(stack));
				if(newstack == NULL){
					printf("Failed to allocate stack in dualtreeAC...WRONG RESULTS\n");
					return;
				}
				newstack->n1 = node1->rptr;
				newstack->n2 = node2;
				newstack->Fbin = Bin2;
				newstack->Lbin = BinStop;
				newstack->next = mystack->next;
				mystack->n1 = node1->lptr;
				mystack->n2 = node2;
				mystack->Fbin = Bin2;
				mystack->Lbin = BinStop;
				mystack->next = newstack;
			}
			continue;
		}

		/* csep is cos(sep) between normal vectors */
		csep = node1->x*node2->x + node1->y*node2->y + node1->z*node2->z;

		/* crange is cos(ang. radius of node1 + ang. radius of node2) */
		crange = node1->clowbound*node2->clowbound - node1->slowbound*node2->slowbound;

		/* Special case if nodes overlap */
		if(crange < csep){
			if(node1->lptr == NULL && node2->lptr == NULL){		/* Both nodes are leafs, time for brute force */
				for(i=BinStart;i<=BinStop;i++){		/* Save the sample counts */
					CntSave[i] = bins[i].Cnt[0];
				}
				
				/* More efficient to avoid references to a node in the inner loop */
				N2Start = node2->Start;
				N2End = node2->End;
				for(i=node1->Start;i<=node1->End;i++){
					for(j=N2Start;j<=N2End;j++){
						csep = DOTP(data[i],data[j]);
						Bin1 = BinStart;
						while(csep > bins[Bin1].limit){
							Bin1++;
						}
						bins[Bin1].Cnt[0]++;
					}
				}
				
				/* Update sample counts */
				Sample1 = node1->Sample;
				Sample2 = node2->Sample;
				if(Sample1 == Sample2){
					for(i=BinStart;i<=BinStop;i++){
						bins[i].Cnt[Sample1] += bins[i].Cnt[0] - CntSave[i];
					}
				}else{
					for(i=BinStart;i<=BinStop;i++){
						Add = bins[i].Cnt[0] - CntSave[i];
						bins[i].Cnt[Sample1] += Add;
						bins[i].Cnt[Sample2] += Add;
					}
				}
				
				/* move on to the next set of nodes in the stack */
				newstack = mystack->next;
				free(mystack);
				mystack = newstack;
				
			}else if(node1->Cnt > node2->Cnt){	/* Explore children of node1 */
				newstack = malloc(sizeof(stack));
				if(newstack == NULL){
					printf("Failed to allocate stack in dualtreeAC...WRONG RESULTS\n");
					return;
				}
				newstack->n1 = node1->rptr;
				newstack->n2 = node2;
				newstack->Fbin = BinStart;
				newstack->Lbin = BinStop;
				newstack->next = mystack->next;
				mystack->n1 = node1->lptr;
				mystack->n2 = node2;
				mystack->Fbin = BinStart;
				mystack->Lbin = BinStop;
				mystack->next = newstack;

			}else{	/* Explore children of node2 */
				newstack = malloc(sizeof(stack));
				if(newstack == NULL){
					printf("Failed to allocate stack in dualtreeAC...WRONG RESULTS\n");
					return;
				}
				newstack->n1 = node1;
				newstack->n2 = node2->rptr;
				newstack->Fbin = BinStart;
				newstack->Lbin = BinStop;
				newstack->next = mystack->next;
				mystack->n1 = node1;
				mystack->n2 = node2->lptr;
				mystack->Fbin = BinStart;
				mystack->Lbin = BinStop;
				mystack->next = newstack;
			}
			continue;
		}

		/* sine of angle between normals */
		ssep = sqrt(1.-csep*csep);
		srange = node1->slowbound*node2->clowbound + node1->clowbound*node2->slowbound;
	
		/* term1 and term2 correspond to those in Dolence & Brunner (2007) */
		term1 = csep*crange;
		term2 = ssep*srange;
		
		mindist = term1 + term2;
		
		/* check if nodes are too distant to be of interest */
		if(mindist < bins[0].limit){
			newstack = mystack->next;
			free(mystack);
			mystack = newstack;
			continue;
		}

		/* Check if nodes include regions 180 deg. apart, otherwise compute maximum separation */
		if(ssep*crange+csep*srange < 0){
			maxdist = -1.;
		}else{
			maxdist = term1-term2;
		}

		/* Bin2 is distance bin for largest separation */
		Bin2 = BinStart;
		while(maxdist > bins[Bin2].limit){
			Bin2++;
		}
		
		/* Bin1 is distance bin for smallest separation */
		Bin1 = Bin2;
		while(mindist > bins[Bin1].limit){
			Bin1++;
		}

		if(Bin1 == Bin2){	/* Bins are the same, so add node1.Cnt*node2.Cnt to total counts and sample counts */
			Add = node1->Cnt*node2->Cnt;
			bins[Bin1].Cnt[0] += Add;
			Sample1 = node1->Sample;
			Sample2 = node2->Sample;
			if(Sample1 == Sample2){
				bins[Bin1].Cnt[Sample1] += Add;
			}else{
				bins[Bin1].Cnt[Sample1] += Add;
				bins[Bin1].Cnt[Sample2] += Add;
			}
		
			/* Move on to the next set of nodes in the stack */
			newstack = mystack->next;
			free(mystack);
			mystack = newstack;
		
		}else if(node1->lptr == NULL && node2->lptr == NULL){	/* nodes are leafs, try single tree subsumption, then brute force */

			BinStart = Bin2;	/* Most distant bin possible */
			BinEnd = Bin1;		/* Closest bin possible */
			
			for(i=BinStart;i<=BinEnd;i++){	/* Save sample counts */
				CntSave[i] = bins[i].Cnt[0];
			}

			N2Start = node2->Start;
			N2End = node2->End;
			for(i=node1->Start;i<=node1->End;i++){
			
				/* do point-node calculations as described in Dolence & Brunner (2007) */
				csep = data[i].x*node2->x + data[i].y*node2->y + data[i].z*node2->z;
				term1 = csep*node2->clowbound;
				term2 = sqrt(1.-csep*csep)*node2->slowbound;
				mindist = term1 + term2;
				maxdist = term1 - term2;
				Bin2 = BinStart;
				while(maxdist > bins[Bin2].limit){
					Bin2++;
				}
				Bin1 = Bin2;
				while(mindist > bins[Bin1].limit){
					Bin1++;
				}
				
				/* Bins are the same, so avoid looking at the individual points in node2 by	*
				 * single-tree subsumption							*/
				if(Bin1 == Bin2){
					bins[Bin1].Cnt[0] += node2->Cnt;
					continue;
				}
				
				/* Brute force calculation necessary */
				/* More efficient to avoid references to a node in the inner loop */
				for(j=N2Start;j<=N2End;j++){
					csep = DOTP(data[i],data[j]);
					Bin1 = Bin2;
					while(csep > bins[Bin1].limit){
						Bin1++;
					}
					bins[Bin1].Cnt[0]++;
				}
			}
			
			/* Update sample counts */
			Sample1 = node1->Sample;
			Sample2 = node2->Sample;
			if(Sample1 == Sample2){
				for(i=BinStart;i<=BinEnd;i++){
					bins[i].Cnt[Sample1] += bins[i].Cnt[0] - CntSave[i];
				}
			}else{
				for(i=BinStart;i<=BinEnd;i++){
					Add = bins[i].Cnt[0] - CntSave[i];
					bins[i].Cnt[Sample1] += Add;
					bins[i].Cnt[Sample2] += Add;
				}
			}
			
			/* Move on to the next set of nodes in the stack */
			newstack = mystack->next;
			free(mystack);
			mystack = newstack;
		
		}else{	/* Explore children */
		
			if(node1->Cnt > node2->Cnt){	/* open node1 */
				newstack = malloc(sizeof(stack));
				if(newstack == NULL){
					printf("Failed to allocate stack in dualtreeAC...WRONG RESULTS\n");
					return;
				}
				newstack->n1 = node1->rptr;
				newstack->n2 = node2;
				newstack->Fbin = Bin2;
				newstack->Lbin = Bin1;
				newstack->next = mystack->next;
				mystack->n1 = node1->lptr;
				mystack->n2 = node2;
				mystack->Fbin = Bin2;
				mystack->Lbin = Bin1;
				mystack->next = newstack;

			}else{	/* open node2 */
				newstack = malloc(sizeof(stack));
				if(newstack == NULL){
					printf("Failed to allocate stack in dualtreeAC...WRONG RESULTS\n");
					return;
				}
				newstack->n1 = node1;
				newstack->n2 = node2->rptr;
				newstack->Fbin = Bin2;
				newstack->Lbin = Bin1;
				newstack->next = mystack->next;
				mystack->n1 = node1;
				mystack->n2 = node2->lptr;
				mystack->Fbin = Bin2;
				mystack->Lbin = Bin1;
				mystack->next = newstack;
			}
		}
	}

	free(CntSave);

	/* All done! */
	
	return;

}

/****************************************************************************************
 *																													 *
 *				      dualtreeCCangular																		 *
 *																													 *
 ****************************************************************************************/
static void dualtreeCCangular(const pcsource data1[],const pcsource data2[],void *nodeA,void *nodeB,bin bins[]){

/* Main routine to cross correlate two data sets A and B.	*
 * Described in Dolence & Brunner (2007).			*/

	unsigned int i,j,N2Start,N2End,Sample1,Sample2,BinEnd;
	unsigned long long int Add,*CntSave;
	double mindist,maxdist,csep,ssep,crange,srange,term1,term2;
	int Bin1,Bin2,BinStart,BinStop;
	angTreeNode *node1,*node2;
	typedef struct stack_{
		angTreeNode *n1,*n2;
		int Fbin,Lbin;
		struct stack_ *next;
	}stack;
	stack *mystack,*newstack;

	BinStart = 0;
	BinStop = NumBins;

	node1 = (angTreeNode *)nodeA;
	node2 = (angTreeNode *)nodeB;
	
	CntSave = malloc((BinStop+2)*sizeof(unsigned long long int));
	if(CntSave == NULL){
		printf("Failed to allocate CntSave in dualtreeCC...WRONG RESULTS\n");
		return;
	}
	
	mystack = malloc(sizeof(stack));
	if(mystack == NULL){
		printf("Failed to allocate stack in dualtreeCC...WRONG RESULTS\n");
		return;
	}
	mystack->n1 = node1;
	mystack->n2 = node2;
	mystack->Fbin = BinStart;
	mystack->Lbin = BinStop;
	mystack->next = NULL;
	
	while(mystack != NULL){	/* Loop until all pairs have been accounted for */
	
		node1 = mystack->n1;
		node2 = mystack->n2;
		BinStart = mystack->Fbin;
		BinStop = mystack->Lbin;

		/* csep is cos(sep) between normal vectors */
		csep = node1->x*node2->x + node1->y*node2->y + node1->z*node2->z;
		
		/* crange is cos(ang. radius of node1 + ang. radius of node2) */
		crange = node1->clowbound*node2->clowbound - node1->slowbound*node2->slowbound;

		/* Special case if nodes overlap */
		if(crange < csep){
			if(node1->lptr == NULL && node2->lptr == NULL){		/* Both nodes are leafs, time for brute force */
				for(i=BinStart;i<=BinStop;i++){		/* Save the sample counts */
					CntSave[i] = bins[i].Cnt[0];
				}
				
				/* More efficient to avoid references to a node in the inner loop */
				N2Start = node2->Start;
				N2End = node2->End;
				for(i=node1->Start;i<=node1->End;i++){
					for(j=N2Start;j<=N2End;j++){
						csep = DOTP(data1[i],data2[j]);
						Bin1 = BinStart;
						while(csep > bins[Bin1].limit){
							Bin1++;
						}
						bins[Bin1].Cnt[0]++;
					}
				}
				
				/* Update sample counts */
				Sample1 = node1->Sample;
				Sample2 = node2->Sample;
				if(Sample1 == Sample2){
					for(i=BinStart;i<=BinStop;i++){
						bins[i].Cnt[Sample1] += bins[i].Cnt[0] - CntSave[i];
					}
				}else{
					for(i=BinStart;i<=BinStop;i++){
						Add = bins[i].Cnt[0] - CntSave[i];
						bins[i].Cnt[Sample1] += Add;
						bins[i].Cnt[Sample2] += Add;
					}
				}
				
				/* move on to the next set of nodes in the stack */
				newstack = mystack->next;
				free(mystack);
				mystack = newstack;
				
			}else if(node1->Cnt > node2->Cnt){	/* Explore children of node1 */
				newstack = malloc(sizeof(stack));
				if(newstack == NULL){
					printf("Failed to allocate stack in dualtreeCC...WRONG RESULTS\n");
					return;
				}
				newstack->n1 = node1->rptr;
				newstack->n2 = node2;
				newstack->Fbin = BinStart;
				newstack->Lbin = BinStop;
				newstack->next = mystack->next;
				mystack->n1 = node1->lptr;
				mystack->n2 = node2;
				mystack->Fbin = BinStart;
				mystack->Lbin = BinStop;
				mystack->next = newstack;

			}else{	/* Explore children of node2 */
				newstack = malloc(sizeof(stack));
				if(newstack == NULL){
					printf("Failed to allocate stack in dualtreeCC...WRONG RESULTS\n");
					return;
				}
				newstack->n1 = node1;
				newstack->n2 = node2->rptr;
				newstack->Fbin = BinStart;
				newstack->Lbin = BinStop;
				newstack->next = mystack->next;
				mystack->n1 = node1;
				mystack->n2 = node2->lptr;
				mystack->Fbin = BinStart;
				mystack->Lbin = BinStop;
				mystack->next = newstack;
			}
			continue;
		}

		ssep = sqrt(1.-csep*csep);
		srange = node1->slowbound*node2->clowbound + node1->clowbound*node2->slowbound;
	
		/* term1 and term2 correspond to those in Dolence & Brunner (2007) */
		term1 = csep*crange;
		term2 = ssep*srange;
		
		/* mindist is cos(minimum separation) */
		mindist = term1 + term2;
		
		/* check if nodes are too distant to be of interest */
		if(mindist < bins[0].limit){
			newstack = mystack->next;
			free(mystack);
			mystack = newstack;
			continue;
		}

		/* Bin1 is distance bin for smallest separation */
		Bin1 = BinStart;
		while(mindist > bins[Bin1].limit){
			Bin1++;
		}

		/* Check if nodes include regions 180 deg. apart, otherwise compute cos(maximum separation) */
		if(ssep*crange+csep*srange < 0){
			maxdist = -1.;
		}else{
			maxdist = term1-term2;
		}

		/* Bin2 is distance bin for largest separation */
		Bin2 = BinStart;
		while(maxdist > bins[Bin2].limit){
			Bin2++;
		}


		if(Bin1 == Bin2){	/* Bins are the same, so add node1.Cnt*node2.Cnt to total counts and sample counts */
			Add = node1->Cnt*node2->Cnt;
			bins[Bin1].Cnt[0] += Add;
			Sample1 = node1->Sample;
			Sample2 = node2->Sample;
			if(Sample1 == Sample2){
				bins[Bin1].Cnt[Sample1] += Add;
			}else{
				bins[Bin1].Cnt[Sample1] += Add;
				bins[Bin1].Cnt[Sample2] += Add;
			}
		
			/* Move on to the next set of nodes in the stack */
			newstack = mystack->next;
			free(mystack);
			mystack = newstack;
		
		}else if(node1->lptr == NULL && node2->lptr == NULL){	/* nodes are leafs, try single tree subsumption, then brute force */

			BinStart = Bin2;	/* Most distant bin possible */
			BinEnd = Bin1;		/* Closest bin possible */
			
			for(i=BinStart;i<=BinEnd;i++){	/* Save sample counts */
				CntSave[i] = bins[i].Cnt[0];
			}
			
			N2Start = node2->Start;
			N2End = node2->End;
			for(i=node1->Start;i<=node1->End;i++){
			
				/* do point-node calculations as described in Dolence & Brunner (2007) */
				csep = data1[i].x*node2->x + data1[i].y*node2->y + data1[i].z*node2->z;
				term1 = csep*node2->clowbound;
				term2 = sqrt(1.-csep*csep)*node2->slowbound;
				mindist = term1 + term2;
				maxdist = term1 - term2;
				Bin1 = BinStart;
				while(mindist > bins[Bin1].limit){
					Bin1++;
				}
				Bin2 = BinStart;
				while(maxdist > bins[Bin2].limit){
					Bin2++;
				}

				/* Bins are the same, so avoid looking at the individual points in node2 by	*
				 * single-tree subsumption							*/
				if(Bin1 == Bin2){
					bins[Bin1].Cnt[0] += node2->Cnt;
					continue;
				}
				
				/* Brute force calculation necessary */
				/* More efficient to avoid references to a node in the inner loop */
				for(j=N2Start;j<=N2End;j++){
					csep = DOTP(data1[i],data2[j]);
					Bin1 = Bin2;
					while(csep > bins[Bin1].limit){
						Bin1++;
					}
					bins[Bin1].Cnt[0]++;
				};
			}
			
			/* Update sample counts */
			Sample1 = node1->Sample;
			Sample2 = node2->Sample;
			if(Sample1 == Sample2){
				for(i=BinStart;i<=BinEnd;i++){
					bins[i].Cnt[Sample1] += bins[i].Cnt[0] - CntSave[i];
				}
			}else{
				for(i=BinStart;i<=BinEnd;i++){
					Add = bins[i].Cnt[0] - CntSave[i];
					bins[i].Cnt[Sample1] += Add;
					bins[i].Cnt[Sample2] += Add;
				}
			}
			
			/* Move on to the next set of nodes in the stack */
			newstack = mystack->next;
			free(mystack);
			mystack = newstack;
		
		}else{	/* Explore children */
		
			if(node1->Cnt > node2->Cnt){	/* open node1 */
				newstack = malloc(sizeof(stack));
				if(newstack == NULL){
					printf("Failed to allocate stack in dualtreeCC...WRONG RESULTS\n");
					return;
				}
				newstack->n1 = node1->rptr;
				newstack->n2 = node2;
				newstack->Fbin = Bin2;
				newstack->Lbin = Bin1;
				newstack->next = mystack->next;
				mystack->n1 = node1->lptr;
				mystack->n2 = node2;
				mystack->Fbin = Bin2;
				mystack->Lbin = Bin1;
				mystack->next = newstack;

			}else{	/* open node2 */
				newstack = malloc(sizeof(stack));
				if(newstack == NULL){
					printf("Failed to allocate stack in dualtreeCC...WRONG RESULTS\n");
					return;
				}
				newstack->n1 = node1;
				newstack->n2 = node2->rptr;
				newstack->Fbin = Bin2;
				newstack->Lbin = Bin1;
				newstack->next = mystack->next;
				mystack->n1 = node1;
				mystack->n2 = node2->lptr;
				mystack->Fbin = Bin2;
				mystack->Lbin = Bin1;
				mystack->next = newstack;
			}
		}
	}

	free(CntSave);

	/* All done! */
	
	return;

}

/* Here to end of file are drivers that call the above functions to get counts. */

#define DTAC(data,node1,node2,bins) (AngOrSpa == 0) ?				\
		dualtreeACangular(data,(node1),(node2),bins) :			\
		dualtreeACspatial(data,(node1),(node2),bins)
		
#define DTCC(data,rand,node1,node2,bins) (AngOrSpa == 0) ?			\
		dualtreeCCangular(data,rand,node1,node2,bins) :			\
		dualtreeCCspatial(data,rand,node1,node2,bins)
		
#define DTAC_ROOT(data,node,child1,child2,root,bins) (AngOrSpa == 0) ?					\
		dualtreeACangular(data,(void *)((((angTreeNode *)node)->child1)->child2),root,bins) :	\
		dualtreeACspatial(data,(void *)((((spaTreeNode *)node)->child1)->child2),root,bins)
		
#define DTCC_ROOT(data,rand,node,child1,child2,root,bins) (AngOrSpa == 0) ?					\
		dualtreeCCangular(data,rand,(void *)((((angTreeNode *)node)->child1)->child2),root,bins) :	\
		dualtreeCCspatial(data,rand,(void *)((((spaTreeNode *)node)->child1)->child2),root,bins)
/*
#define DTAC(data,child1,child2,root,bins) ((*dualtreeAC)==dualtreeACangular) ?			\
		(*dualtreeAC)(data,(void *)(atemp->child1)->child2,root,bins) :			\
		(*dualtreeAC)(data,(void *)(stemp->child1)->child2,root,bins)

#define DTCC(data,rand,child1,child2,root,bins) ((*dualtreeCC)==dualtreeCCangular) ?		\
		(*dualtreeCC)(data,rand,(void *)(atemp->child1)->child2,root,bins) :		\
		(*dualtreeCC)(data,rand,(void *)(stemp->child1)->child2,root,bins)
*/
#ifndef USE_MPI

#ifndef USE_OMP
void ac_serial(const pcsource data[],void *root1,void *root2,bin bins[]){

	TIMESTART(startTime);
	DTAC(data,root1,root2,bins);
	TIMESTOP(wallTimeAC,startTime);
	return;
}

void cc_serial(const pcsource dataA[],const pcsource dataB[],void *root1,void *root2,bin bins[]){

	TIMESTART(startTime);
	DTCC(dataA,dataB,root1,root2,bins);
	TIMESTOP(wallTimeCC,startTime);
	return;
}
#endif

#ifdef USE_OMP
static int GlobalPos;

void ac_shared(const pcsource data[],void *worknodes[],void *root,bin bins[],int NumWorkNodes,int MyId){

	int MyNode;
	clock_t stopTime,diffTime;

	TIMESTART(startTime);
	
	MyNode = MyId;

	#pragma omp master
	{
		GlobalPos = NumThreads;
	}
	#pragma omp barrier

	while(MyNode < NumWorkNodes){
		DTAC(data,worknodes[MyNode],root,bins);
		#pragma omp critical                              
		{
			MyNode = GlobalPos;
			GlobalPos++;
		}
	}

	#pragma omp barrier
	
	TIMESTART(stopTime);
	#ifdef TIMING
	diffTime = difftime(stopTime,startTime);
	#pragma omp atomic
	wallTimeAC += diffTime;
	#endif

	return;
}
	
void cc_shared(const pcsource dataA[],const pcsource dataB[],void *worknodes[],void *root,bin bins[],int NumWorkNodes,int MyId){

	int MyNode;
	clock_t stopTime,diffTime;

	TIMESTART(startTime);
	
	MyNode = MyId;

	#pragma omp master
	{
		GlobalPos = NumThreads;
	}
	#pragma omp barrier

	while(MyNode < NumWorkNodes){
		DTCC(dataA,dataB,worknodes[MyNode],root,bins);
		#pragma omp critical
		{
			MyNode = GlobalPos;
			GlobalPos++;
		}
	}

	#pragma omp barrier
	
	TIMESTART(stopTime);
	#ifdef TIMING
	diffTime = difftime(stopTime,startTime);
	#pragma omp atomic
	wallTimeCC += diffTime;
	#endif

	return;
}
#endif
#endif

#ifdef USE_MPI
static int GlobalPos,RootPos,NumTerminated;

static void processRequests(int TestPos){

/* process all messages, responding with the next worknode or the terminate message */

	int IsMessage,Source,Terminate;
	MPI_Status status;
	#ifdef TIMING
	double begin;
	#endif
	
	TIMECHECK(begin);

	Terminate = -1;

	MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&IsMessage,&status);
	while(IsMessage != FALSE){      /* Process messages until queue is empty */
		MPI_Recv(&Source,1,MPI_INT,status.MPI_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
		#pragma omp critical
		{
			if(GlobalPos < TestPos){  /* Send the next node in the list */
				MPI_Send(&GlobalPos,1,MPI_INT,status.MPI_SOURCE,Source,MPI_COMM_WORLD);
				GlobalPos++;
			}else{  /* Out of work, send a terminate message */
				MPI_Send(&Terminate,1,MPI_INT,status.MPI_SOURCE,Source,MPI_COMM_WORLD);
				NumTerminated++;
			}
		}
		MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&IsMessage,&status);       /* Check for messages*/
	}

	TIMEDIFF(wallComm,begin);
	
	return;
}
#endif


#if defined (USE_MPI) && defined (USE_OMP)

void ac_dist_shared(const pcsource data[],void *worknodes[],void *root,bin bins[],int NumWorkNodes,int MyId){


	int MyNode[2],MyNodeNow,Terminate,Num2Term,MyProc,Flag;
	double stopTime,diffTime;
	MPI_Status statuses[2];
	MPI_Request request[2];

	TIMECHECK(startTime);
	
	DEBUGPRINT("Entering ac_dist_shared\n");
	
	if(NumProcs == 1 && NumThreads == 1){
		DTAC(data,root,root,bins);
	}else if(NumProcs == 1 && NumThreads != 1){
		GlobalPos = NumThreads;
		MyNodeNow = MyId;
		while(MyNodeNow < NumWorkNodes){
			DTAC(data,worknodes[MyNodeNow],root,bins);
			#pragma omp critical
			{
				MyNodeNow = GlobalPos;
				GlobalPos++;
			}
		}
	}else{

	Num2Term = NumThreads*(NumProcs - 1);
	MyProc = MyRank*NumThreads + MyId;
	Terminate = -1;

	if(MyRank == 0){
		if(MyId == 0){
			#pragma omp critical
			{
				RootPos = NumWorkNodes-1;
				GlobalPos = 2*Num2Term+NumThreads-1;
			}
			NumTerminated = 0;
			while(NumTerminated != Num2Term){
				processRequests(RootPos+1);
				if(GlobalPos <= RootPos){
					DTAC_ROOT(data,worknodes[RootPos],lptr,lptr,root,bins);
					processRequests(RootPos);
					DTAC_ROOT(data,worknodes[RootPos],lptr,rptr,root,bins);
					processRequests(RootPos);
					DTAC_ROOT(data,worknodes[RootPos],rptr,lptr,root,bins);
					processRequests(RootPos);
					DTAC_ROOT(data,worknodes[RootPos],rptr,rptr,root,bins);
					#pragma omp critical
					{ RootPos--; }
				}
			}
		}else{
			Flag = 1;
			MyNodeNow = MyId-1;
			while(Flag == 1){
				DTAC(data,worknodes[MyNodeNow],root,bins);
				#pragma omp critical
				{
					MyNodeNow = GlobalPos;
					GlobalPos++;
					if(MyNodeNow >= RootPos) Flag = 0;
				}
			}
		}
	}else{
		/* Figure out first two nodes to work with */
		MyNode[0] = 2*MyProc - NumThreads - 1;
		MyNode[1] = MyNode[0]+1;

		while(MyNode[1] != Terminate){ /* Loop to work on my nodes and get new nodes as needed */

			/* Ask for the next node */
			#pragma omp critical
			{
			MPI_Isend(&MyProc,1,MPI_INT,0,MyProc,MPI_COMM_WORLD,&request[0]);
			}
			MyNodeNow = MyNode[0];
			MyNode[0] = MyNode[1];

			/* Receive new node that I am responsible for */
			#pragma omp critical
			{
			MPI_Irecv(&MyNode[1],1,MPI_INT,0,MyProc,MPI_COMM_WORLD,&request[1]);
			}
			/* Do work on the first node in my queue */
			DTAC(data,worknodes[MyNodeNow],root,bins);

			/* Wait for all communication to finish */
			#pragma omp critical
			{
			MPI_Waitall(2,request,statuses);
			}

		}
		/* Do the last one in my queue */
		DTAC(data,worknodes[MyNode[0]],root,bins);
	}
	}

	#pragma omp barrier
	DEBUGPRINT("Exiting ac_dist_shared\n");

	TIMECHECK(stopTime);
	#ifdef TIMING
	diffTime = stopTime-startTime;
	#pragma omp atomic
	wallTimeAC += diffTime;
	#endif
	
	return;
}

void cc_dist_shared(const pcsource data[],const pcsource rand[],void *worknodes[],void *root,bin bins[],int NumWorkNodes,int MyId){

/* Routine to handle the parallel distribution of work in a cross correlation.	*
 * Described in Dolence & Brunner (2007).					*/

	int MyNode[2],MyNodeNow,Terminate,Num2Term,MyProc;
	double stopTime,diffTime;
	MPI_Status statuses[2];
	MPI_Request request[2];

	TIMECHECK(startTime);
	
	DEBUGPRINT("Entering cc_dist_shared\n");

        if(NumProcs == 1 && NumThreads == 1){
                for(MyNodeNow=0;MyNodeNow<NumWorkNodes;MyNodeNow++){
			DTCC(data,rand,worknodes[MyNodeNow],root,bins);
		}
        }else if(NumProcs == 1 && NumThreads != 1){
                GlobalPos = NumThreads;
                MyNodeNow = MyId;
                while(MyNodeNow < NumWorkNodes){
                        DTCC(data,rand,worknodes[MyNodeNow],root,bins);
                        #pragma omp critical
                        {
                                MyNodeNow = GlobalPos;
                                GlobalPos++;
                        }
                }
        }else{
	
	Terminate = -1;
	Num2Term = NumThreads*(NumProcs - 1);
	MyProc = MyRank*NumThreads + MyId;
	

	if(MyRank == 0){
		#pragma omp master
		{ GlobalPos = 2*Num2Term + NumThreads; }	/* Next node that hasn't been assigned to a process */
		#pragma omp barrier
		if(MyId == 0){
			NumTerminated = 0;
			RootPos = 0;				/* My first node to work on */
			while(NumTerminated != Num2Term){	/* Loop until everyone knows to quit */
				processRequests(NumWorkNodes);
				if(RootPos < NumWorkNodes){	/* if I have a node to work on, do so now */
					DTCC_ROOT(data,rand,worknodes[RootPos],lptr,lptr,root,bins);
					processRequests(NumWorkNodes);
					DTCC_ROOT(data,rand,worknodes[RootPos],lptr,rptr,root,bins);
					processRequests(NumWorkNodes);;
					DTCC_ROOT(data,rand,worknodes[RootPos],rptr,lptr,root,bins);
					processRequests(NumWorkNodes);
					DTCC_ROOT(data,rand,worknodes[RootPos],rptr,rptr,root,bins);
					#pragma omp critical
					{
						RootPos = GlobalPos;
						GlobalPos++;
					}
				}
			}
			/* Terminate messages have been sent to everyone */
			if(RootPos < NumWorkNodes){	/* Do work if I've got some left over */
				DTCC(data,rand,worknodes[RootPos],root,bins);
			}
		}else{
			MyNodeNow = MyId;
			while(MyNodeNow < NumWorkNodes){
				DTCC(data,rand,worknodes[MyNodeNow],root,bins);
				#pragma omp critical
				{
					MyNodeNow = GlobalPos;
					GlobalPos++;
				}
			}
		}	
	}else{
		/* Figure out first two nodes to work on */
		MyNode[0] = 2*MyProc-NumThreads;
		MyNode[1] = MyNode[0]+1;
		
		while(MyNode[1] != Terminate){	/* Loop until there is no more work to do */
			#pragma omp critical
			{ MPI_Isend(&MyProc,1,MPI_INT,0,MyProc,MPI_COMM_WORLD,&request[0]); }
			MyNodeNow = MyNode[0];
			MyNode[0] = MyNode[1];
			#pragma omp critical
			{ MPI_Irecv(&MyNode[1],1,MPI_INT,0,MyProc,MPI_COMM_WORLD,&request[1]); }
			DTCC(data,rand,worknodes[MyNodeNow],root,bins);
			#pragma omp critical
			{ MPI_Waitall(2,request,statuses); }
		}
		/* Do the last one in my queue */
		DTCC(data,rand,worknodes[MyNode[0]],root,bins);
	}
	}

	#pragma omp barrier

	DEBUGPRINT("Exiting cc_dist_shared\n");

	TIMECHECK(stopTime);
	#ifdef TIMING
	diffTime = stopTime-startTime;
	#pragma omp atomic
	wallTimeCC += diffTime;
	#endif
	
	return;
}
#endif


#ifdef USE_MPI
#ifndef USE_OMP

void ac_dist(const pcsource data[],void *worknodes[],void *root,bin bins[],int NumWorkNodes){

/* Routine to handle the parallel distribution of work in an autocorrelation.	*
 * Described in Dolence & Brunner (2007).					*/

	int MyNodeNow,MyNode[2],Terminate,NPm1;
	MPI_Status statuses[2];
	MPI_Request request[2];

	TIMESTART(startTime);
	
	DEBUGPRINT("Entering ac_dist\n");

	Terminate = -1;
	NPm1 = NumProcs-1;

	/* No need to break up the problem if only one process is available */
	if(NumProcs == 1){
		DTAC(data,root,root,bins);
	}else{
	
	if(MyRank == 0){
		NumTerminated = 0;
		GlobalPos = 2*NPm1;	/* Next node to send to another process */
		RootPos = NumWorkNodes-1;	/* My first node, starting at the end of the list */
		while(NumTerminated != NPm1){	/* Loop until everyone knows to quit */
			processRequests(RootPos+1);
			if(GlobalPos <= RootPos){
				DTAC_ROOT(data,worknodes[RootPos],lptr,lptr,root,bins);
				processRequests(RootPos);
				DTAC_ROOT(data,worknodes[RootPos],lptr,rptr,root,bins);
				processRequests(RootPos);
				DTAC_ROOT(data,worknodes[RootPos],rptr,lptr,root,bins);
				processRequests(RootPos);
				DTAC_ROOT(data,worknodes[RootPos],rptr,rptr,root,bins);
				RootPos--;
			}
		}	
	}else{

		/* Figure out first two nodes to work with */
		MyNode[0] = 2*(MyRank-1);
		MyNode[1] = MyNode[0]+1;
		
		while(MyNode[1] != Terminate){ /* Loop to work on my nodes and get new nodes as needed */

			/* Ask for the next node */
			MPI_Isend(&MyRank,1,MPI_INT,0,MyRank,MPI_COMM_WORLD,&request[0]);
			MyNodeNow = MyNode[0];
			MyNode[0] = MyNode[1];

			/* Receive new node that I am responsible for */
			MPI_Irecv(&MyNode[1],1,MPI_INT,0,MyRank,MPI_COMM_WORLD,&request[1]);

			/* Do work on the first node in my queue */
			DTAC(data,worknodes[MyNodeNow],root,bins);

			/* Wait for all communication to finish */
			MPI_Waitall(2,request,statuses);
		}
		/* Do the last one in my queue */
		DTAC(data,worknodes[MyNode[0]],root,bins);
	}
	}

	DEBUGPRINT("Exiting ac_dist\n");

	TIMESTOP(wallTimeAC,startTime);
	
	return;
}

void cc_dist(const pcsource data[],const pcsource rand[],void *worknodes[],void *root,bin bins[],int NumWorkNodes){

/* Routine to handle the parallel distribution of work in a cross correlation.	*
 * Described in Dolence & Brunner (2007).					*/

	int MyNode[2],MyNodeNow,Terminate,NPm1;
	MPI_Status statuses[2];
	MPI_Request request[2];

	TIMESTART(startTime);
	
	DEBUGPRINT("Entering cc_dist\n");

	Terminate = -1;
	NPm1 = NumProcs-1;

	if(NumProcs == 1){
		for(GlobalPos = 0;GlobalPos<NumWorkNodes;GlobalPos++){
			DTCC(data,rand,worknodes[GlobalPos],root,bins);
		}
	}else{

	if(MyRank == 0){
		NumTerminated = 0;
		GlobalPos = 2*NumProcs-1;	/* Next node that hasn't been assigned to a process */
		RootPos = 0;	/* My first node to work on */
		while(NumTerminated != NPm1){	/* Loop until everyone knows to quit */
			processRequests(NumWorkNodes);
			if(RootPos < NumWorkNodes){	/* if I have a node to work on, do so now */
				DTCC_ROOT(data,rand,worknodes[RootPos],lptr,lptr,root,bins);
				processRequests(NumWorkNodes);
				DTCC_ROOT(data,rand,worknodes[RootPos],lptr,rptr,root,bins);
				processRequests(NumWorkNodes);
				DTCC_ROOT(data,rand,worknodes[RootPos],rptr,lptr,root,bins);;
				processRequests(NumWorkNodes);
				DTCC_ROOT(data,rand,worknodes[RootPos],rptr,rptr,root,bins);
				RootPos = GlobalPos;
				GlobalPos++;
			}
		}
		/* Terminate messages have been sent to everyone */
		if(RootPos < NumWorkNodes){	/* Do work if I've got some left over */
			DTCC(data,rand,worknodes[RootPos],root,bins);
		}	
	}else{
		/* Figure out first two nodes to work on */
		MyNode[0] = 2*MyRank-1;
		MyNode[1] = MyNode[0]+1;
		
		while(MyNode[1] != Terminate){	/* Loop until there is no more work to do */
			MPI_Isend(&MyRank,1,MPI_INT,0,MyRank,MPI_COMM_WORLD,&request[0]);
			MyNodeNow = MyNode[0];
			MyNode[0] = MyNode[1];
			MPI_Irecv(&MyNode[1],1,MPI_INT,0,MyRank,MPI_COMM_WORLD,&request[1]);
			DTCC(data,rand,worknodes[MyNodeNow],root,bins);
			MPI_Waitall(2,request,statuses);
		}
		/* Do the last one in my queue */
		DTCC(data,rand,worknodes[MyNode[0]],root,bins);
	}
	}

	DEBUGPRINT("Exiting cc_dist\n");
	
	TIMESTOP(wallTimeCC,startTime);

	return;
}
#endif
#endif

