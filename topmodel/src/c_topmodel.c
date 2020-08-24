#define	MAIN
#include "topmodel.h"
#undef	MAIN

void c_topmodel(double *parameters,
	     double *topidx,
             double *delay,
	     double *rain,
	     double *ETp,
	     double *Qobs,
	     int *nidxclass,
	     int *ntimestep,
             int *iterations,
             int *nch,
	     int *verbose,
	     double *result)
{
	int i,j;

	topidx_calc(topidx, *nidxclass);
	memory_allocation(*nch, *ntimestep, *nidxclass);
	
	if(*iterations > 1) Rprintf("Iteration:         ");
	
	#ifdef win32
	   R_flushConsole();
	   R_ProcessEvents();
	#endif

	for(i=0;i<*iterations;i++) {

		R_CheckUserInterrupt();
		if(*iterations > 1) Rprintf ("\b\b\b\b\b\b\b\b%8i",i+1);

		param_init(parameters, delay, *nch, i, *nidxclass, *ntimestep);

		/* run the model for each time step */

		for(j=0; j<*ntimestep; j++)
			run_topmodel(rain,ETp,*nidxclass,j,*ntimestep);

		/* write outputs */

		output(Qobs, result, *ntimestep, *iterations, *verbose, *nidxclass, i);
	
	}

	if(*iterations > 1) Rprintf("\n");
	memory_free(*nch, *ntimestep, *nidxclass);
	return;
}
