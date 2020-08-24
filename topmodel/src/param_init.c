#include "topmodel.h"

void param_init(double *parameters, double *delay, int nch, int iteration, int nidxclass, int ntimestep)
{
	int nparam = 11;
	int i;
	double A;

	params.qs0    = parameters[0 + nparam*iteration];
	params.lnTe   = parameters[1 + nparam*iteration];
	params.m      = parameters[2 + nparam*iteration];
	params.Sr0    = parameters[3 + nparam*iteration];
	params.Srmax  = parameters[4 + nparam*iteration];
	params.td     = parameters[5 + nparam*iteration];
	params.vch    = parameters[6 + nparam*iteration];
	params.vr     = parameters[7 + nparam*iteration];
	params.K0     = parameters[8 + nparam*iteration];
	params.CD     = parameters[9 + nparam*iteration];
	params.dt     = parameters[10+ nparam*iteration];

/* reading delay function */

	for(i=0; i<nch; i++){
		params.d[i] = delay[i];
		params.Ad_r[i] = delay[i+nch];
	}

	misc.lambda = get_lambda(nidxclass);
	misc.lnTe   = params.lnTe + log(params.dt);
	misc.vch    = params.vch * params.dt;
	misc.vr     = params.vr * params.dt;
	misc.qs0    = params.qs0 * params.dt;
	misc.qss    = exp(misc.lnTe - misc.lambda);

/* Initialisation functions:
 *
 * 1. LnTe (Areal average of ln(T0) = eq. 6.21c):
 *    ln(Te*dt) = lnTe + log(dt)
 *
 * 2. qss (discharge when S_mean is zero, i.e. saturation:
 *    qss = exp(-gamma)   (see Beven 2000 p. 212)
 *    and gamma = lambda - lnTe
 */

	get_Ad(nch); 

/* initialisation of S_mean by means of eq. 6.33 */

	misc.S_mean[0] = - params.m * log(misc.qs0 / misc.qss);

/* initialisation of Qt 
   Qt = qs0 in the beginning before ndelay
   and has recession curve between ndelay and (ndelay+reach) */

	for(i=0; i<ntimestep; i++) misc.Qt[i] = 0.0;
	for(i=0; i<misc.ndelay; i++) misc.Qt[i] = misc.qs0;

	A = 0.0;

	for(i=0; i<misc.nreach; i++){
		A += misc.Ad[i];
		misc.Qt[misc.ndelay + i] = misc.qs0 * (1 - A);
	}
	
/* initialisation of Srz and Suz for the first timestep */	

	for(i=0; i<nidxclass; i++){
		misc.Srz[0][i] = params.Sr0;
		misc.Suz[0][i] = 0.0;
	}

}
