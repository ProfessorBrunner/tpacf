
/* define some global variables */
int NumBins,NumSamples,AngOrSpa;
#ifdef USE_MPI
int MyRank,NumProcs;
#endif
#ifdef USE_OMP
int NumThreads;
#endif

#ifdef TIMING
double wallTimeAC,wallTimeCC,wallTimeTotal;
double wallDataRead,wallTreeRead,wallSetup;
#ifdef USE_OMP
double wallReduceOmp;
#endif
#ifdef USE_MPI
double startTime,wallReduce,wallComm;
#else
clock_t startTime,endTime;
#endif
#endif
/* end globals */
