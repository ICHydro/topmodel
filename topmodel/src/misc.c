#include "topmodel.h"

double get_lambda(int nidxclass)
{
  int    i;
  double retval;

  /* take average because idx values represent lower limit of the interval! */

  retval = 0.0;
  for(i=1; i<nidxclass; i++) {
    retval += idxstats.Aatb_r[i] * (idxstats.atb[i] + idxstats.atb[i - 1]) / 2.0;
  }
  return retval;
}

void get_Ad(int nch)
{
  int i, j, t;
  double A1, A2;

  misc.tch = (double *) R_alloc(nch, sizeof(double));
  misc.tch[0] = params.d[0] / misc.vch;

  for(i=1; i<nch; i++)
    misc.tch[i] = misc.tch[0] + (params.d[i]-params.d[0]) / misc.vr;

  /* tch = number of timesteps for a each Aatb to reach the outlet (not yet discretized) */

  misc.nreach = (int) misc.tch[nch - 1];
  if((double) misc.nreach < misc.tch[nch - 1]) misc.nreach++;

  /* nreach = number of timesteps for the farthest part of the catchment (discretized and rounded upwards) */

  misc.ndelay = (int) misc.tch[0];

  /* delay from outlet to catchment to gauge, rounded downwards, discrete */

  misc.nreach -= misc.ndelay;

  /* delay within catchment, discretized */

  misc.Ad = (double *) R_alloc(misc.nreach, sizeof(double));

  for(i=0; i<misc.nreach; i++){
    t = misc.ndelay + i + 1;
    if(t > misc.tch[nch - 1]){
      misc.Ad[i] = 1.0;
    }else{
      for(j=1; j<nch; j++){
	if(t <= misc.tch[j]){
	  misc.Ad[i] = params.Ad_r[j - 1] +
	    (params.Ad_r[j] - params.Ad_r[j - 1]) *
	    (t - misc.tch[j - 1]) / (misc.tch[j] - misc.tch[j - 1]);
	  break;
	}
      }
    }
  }

  /* builds a *cumulative* area vector over the total delay of the catchment, to determine
     which timestep gets which part of the outflow
     RESULT = Ad */

  A1 = misc.Ad[0];
  for(i=1; i<misc.nreach; i++){
    A2 = misc.Ad[i];
    misc.Ad[i] = A2 - A1;
    A1 = A2;
  }

/* now Ad not cumulative any more */

}



