#include "topmodel.h"

void output(double *Qobs, double *result, int ntimestep, int iterations, int verbose, int nidxclass, int i)
{
	int start, j;

/* do we return simulated discharge...? */

	if(Qobs[0] == -9999) {				
		for(j=0; j < ntimestep; j++){
			result[i * ntimestep + j] = misc.Qt[j];
		}
		if(verbose == 6) {
			/* output qo */
			start = ntimestep * iterations;
			for(j=0; j < ntimestep; j++){
				result[start + i * ntimestep+j] = misc.qo[j][nidxclass];
			}
			/* output qs */
			start = 2 * ntimestep * iterations;
			for(j=0; j < ntimestep; j++){
				result[start + i * ntimestep+j] = misc.qs[j];
			}
			/* output S_mean */
			start = 3 * ntimestep * iterations;
			for(j=0; j < ntimestep; j++){
				result[start + i * ntimestep+j] = misc.S_mean[j];
			}
			/* output infiltration excess (fex) */
			start = 4 * ntimestep * iterations;
			for(j=0; j < ntimestep; j++){
				result[start + i * ntimestep+j] = misc.fex[j];
			}
			/* output actual ET */
			start = 5 * ntimestep * iterations;
			for(j=0; j < ntimestep; j++){
				result[start + i * ntimestep+j] = misc.Ea[j][nidxclass];
			}
		}
	}

/*  ...or do we calculate directly the Nash Sutcliffe efficiency? */

	else {						
		result[i] = get_Eff(misc.Qt, Qobs, ntimestep);
	}
}

