#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <R.h>

#define	ZERO		0.0000001

/* output.c */
void output(double *, double *, int, int, int, int, int);
/* start.c */
void memory_allocation(int, int, int);
void memory_free(int, int, int);
void topidx_calc(double *, int);
/* NS.c */
double get_Eff(double *, double *, int);
/* misc */
void get_Ad(int);
double get_lambda(int);
/* param_init */
void param_init(double *, double *, int, int, int, int);
/* get_f.c */
double get_f(double, double, double, double, double, double);
/* cores */
void run_topmodel(double *, double *, int, int, int);

#ifdef MAIN
#	define	GLOBAL
#else
#	define	GLOBAL	extern
#endif

/* Topographic index statistics */
GLOBAL	struct
{
	double	*atb, *Aatb_r;
} idxstats;

GLOBAL	struct
{
	double	qs0, lnTe, m, Sr0, Srmax, td, vch, vr, n;
	double	K0, CD, dt;
	/* channel parameters */
	double	*d, *Ad_r;
} params;

GLOBAL	struct
{
	/* Model efficiency */
	double	Em;
	int	ndelay, nreach;
	double	lnTe, vch, vr;
	double	lambda;
	double	qss, qs0;
	double	Qobs_mean;
	/* params.nch's */
	double	*tch;
	/* misc.nreach's */
	double	*Ad;
	/* input.ntimestep's */
	double	*Qobs;
	double	*Qt;
	double	*qs;	/* spatially constant! */
	double	*S_mean;
	double	*f;
	double	*fex;
	/* input.ntimestep * (nidxclass + 1)'s */
	double	**qt, **qo, **qv, **qint;
	/* input.ntimestep * nidxclass's */
	double	**Srz, **Suz;
	double	**S;
	double	**Ea;
	double	**ex;
	/* Miscellaneous variables */
	int	timestep, idxclass;
} misc;



