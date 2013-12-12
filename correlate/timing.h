
#ifndef USE_MPI
	#ifdef TIMING
		#include <time.h>
	#endif
#endif


/* define some timing macros and variable defs*/
#ifdef TIMING
	extern double wallTimeAC,wallTimeCC,wallTimeTotal;
	extern double wallDataRead,wallTreeRead,wallSetup;
	#ifdef USE_OMP
		extern double wallReduceOmp;
	#endif
	#ifdef USE_MPI
		extern double startTime,wallReduce,wallComm;
		#define TIMESTART(startTimeVar) 				\
			MPI_Barrier(MPI_COMM_WORLD);				\
			startTimeVar = MPI_Wtime()		
		#define TIMESTOP(totalTimeVar,startTimeVar)			\
			MPI_Barrier(MPI_COMM_WORLD);				\
			totalTimeVar += MPI_Wtime() - (startTimeVar)
		#define TIMECHECK(timeVar) timeVar = MPI_Wtime()
		#define TIMEDIFF(timeSumVar,timeVar) timeSumVar += MPI_Wtime() - (timeVar)
	#else
		extern clock_t startTime,endTime;
		#define TIMESTART(startTimeVar) time(&(startTimeVar))
		#define TIMESTOP(totalTimeVar,startTimeVar)			\
			time(&endTime);						\
			totalTimeVar += difftime(endTime,(startTimeVar))
	#endif
#else
	#define TIMECHECK(timeVar)
	#define TIMEDIFF(timeSumVar,timeVar)
	#define TIMESTART(startTimeVar)
	#define TIMESTOP(totalTimeVar,startTimeVar)
#endif

void reportTimings(void);

