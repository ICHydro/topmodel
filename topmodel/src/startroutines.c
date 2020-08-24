#include "topmodel.h"

void topidx_calc(double *topidx, int nidxclass)
{
	int	i, j;
	double	x;

	idxstats.atb = (double *) Calloc(nidxclass, double);
	idxstats.Aatb_r = (double *) Calloc(nidxclass, double);

	for(i=0; i<nidxclass; i++) {
		idxstats.atb[i] = topidx[i];
		idxstats.Aatb_r[i] = topidx[i+nidxclass];
	}

/* sort index classes from high to low */

	for(i=0; i<nidxclass; i++){
		for(j=i; j<nidxclass; j++){
			if(idxstats.atb[i] < idxstats.atb[j]){
				x = idxstats.atb[i];
				idxstats.atb[i] = idxstats.atb[j];
				idxstats.atb[j] = x;
				x = idxstats.Aatb_r[i];
				idxstats.Aatb_r[i] = idxstats.Aatb_r[j];
				idxstats.Aatb_r[j] = x;
			}
		}
	}

	return;
}


void memory_allocation(int nch, int ntimestep, int nidxclass)
{
	int i;

	misc.Qt = (double *) Calloc(ntimestep, double);
	misc.S_mean = (double *) Calloc(ntimestep, double);

	params.d = (double *) Calloc(nch, double);
	params.Ad_r = (double *) Calloc(nch, double);

	misc.Srz = (double **) Calloc(ntimestep, double *); /* Root zone storage deficit */
	misc.Suz = (double **) Calloc(ntimestep, double *); /* Unsaturated zone storage */

	misc.S = (double **) Calloc(ntimestep, double *);
	misc.Ea = (double **) Calloc(ntimestep, double *);
	misc.ex = (double **) Calloc(ntimestep, double *);

	misc.qt = (double **) Calloc(ntimestep, double *);
	misc.qo = (double **) Calloc(ntimestep, double *);
	misc.qv = (double **) Calloc(ntimestep, double *);
	misc.qint = (double **) Calloc(ntimestep, double *);

	misc.qs = (double *) Calloc(ntimestep, double);
	misc.f = (double *) Calloc(ntimestep, double);
	misc.fex = (double *) Calloc(ntimestep, double);

	for(i=0; i<ntimestep; i++){
		misc.Srz[i] = (double *) Calloc(nidxclass, double);
		misc.Suz[i] = (double *) Calloc(nidxclass, double);

		misc.S[i]  = (double *) Calloc(nidxclass, double);
		misc.Ea[i] = (double *) Calloc((nidxclass + 1), double);
		misc.ex[i] = (double *) Calloc((nidxclass + 1), double);

		misc.qt[i] = (double *) Calloc((nidxclass + 1), double);
		misc.qo[i] = (double *) Calloc((nidxclass + 1), double);
		misc.qv[i] = (double *) Calloc((nidxclass + 1), double);
		misc.qint[i] = (double *) Calloc((nidxclass + 1), double);
	}

	return;
}

void memory_free(int nch, int ntimestep, int nidxclass)
{
	int i;

	for(i=0; i<ntimestep; i++){
		Free(misc.Srz[i]);
		Free(misc.Suz[i]);

		Free(misc.S[i]);
		Free(misc.Ea[i]);
		Free(misc.ex[i]);

		Free(misc.qt[i]);
		Free(misc.qo[i]);
		Free(misc.qv[i]);
		Free(misc.qint[i]);
	}

	Free(misc.Qt);
	Free(misc.S_mean);

	Free(params.d);
	Free(params.Ad_r);

	Free(misc.Srz);
	Free(misc.Suz);

	Free(misc.S);
	Free(misc.Ea);
	Free(misc.ex);

	Free(misc.qt);
	Free(misc.qo);
	Free(misc.qv);
	Free(misc.qint);

	Free(misc.qs);
	Free(misc.f);
	Free(misc.fex);

	Free(idxstats.atb);
	Free(idxstats.Aatb_r);

/* a few missing! */

	return;
}
