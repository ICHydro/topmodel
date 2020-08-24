#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <R.h>

#define	ZERO		0.0000001

/* output.c */
void	output();
/* start.c */
void	memory_allocation();
void	memory_free();
void	topidx_calc();
/* NS.c */
double	get_Eff();
/* misc */
void	get_Ad();
double	get_lambda();
/* param_init */
void	param_init();
/* infiltration.c */
double	get_f();
/* cores */
void	run_topmodel();
void    run_topsat();

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



