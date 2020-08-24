#include "topmodel.h"

void run_topmodel(double *rain, double *ET0, int nidxclass, int i, int ntimestep)
{
  int    j, k;
  double Aatb_r, _qo, _qv;

  /* initialise the fluxes */

  misc.qt[i][nidxclass] = 0.0;
  misc.qo[i][nidxclass] = 0.0;
  misc.qv[i][nidxclass] = 0.0;
  misc.qs[i] = 0.0;
  misc.f[i] = rain[i];               /* By default all rain infiltrates */
  misc.fex[i] = 0.0;                 /* and therefore fex is zero */

  /* calculate infiltration and redirect any excess infiltration to fex */

  if(params.infex){
    misc.f[i] = params.dt * get_f((i + 1) * params.dt, rain[i] / params.dt,
				  params.CD, params.K0, params.m, params.dt);
    if(misc.f[i]<0) misc.f[i] = rain[i];
    /* necessary? -> yes! but would be good to find out why ...*/
    misc.fex[i] = rain[i] - misc.f[i];
  }

  /* Srz = Root zone storage deficit
     Suz = Unsaturated (gravity drainage) zone storage */

  if(i){
    for(j=0; j<nidxclass; j++){
      misc.Srz[i][j] = misc.Srz[i-1][j];
      misc.Suz[i][j] = misc.Suz[i-1][j];
    }
  }

  misc.qs[i] = misc.qss * exp(- misc.S_mean[i] / params.m);	/* eq. 18.16 */

/* qs = Subsurface flow per unit area
   qss = saturated zone flow = exp(-gamma) = exp(lnTe-lambda) */

/*******************************
 * begin loop on index classes *
 *******************************/

  for(j=0; j<nidxclass; j++){

    /* Area index class = (area this class + area next class)/2 (except for last, there just last area/2)
       This is because the topidx values represent the LOWER limit respective area fraction */

    Aatb_r = (idxstats.Aatb_r[j] + 
	      (j < nidxclass - 1 ? idxstats.Aatb_r[j + 1] : 0.0)) / 2.0;

    /* calculate local storage deficit		
       ! In topmodel T0 is assumed to be spatially constant so (lnT0 - lnTe) = 0        */

    misc.S[i][j] = misc.S_mean[i] + params.m * (misc.lambda - idxstats.atb[j]); /* (eq. 18.8) */

/* rain first enters the root zone (written as deficit)	*/

    if(misc.S[i][j] < 0.0) misc.S[i][j] = 0.0;
    misc.Srz[i][j] -= misc.f[i];

    /* only when root zone is filled, water goes to the unsaturated zone
       rz is deficit but uz not, therefore substraction -> a negative Srz adds to Suz */

    if(misc.Srz[i][j] < 0.0){
      misc.Suz[i][j] -= misc.Srz[i][j];
      misc.Srz[i][j] = 0.0;
    }

/* if Suz exceeds S, then excess flow (ex) occurs (saturated overland flow)
   S determines the storage capacity of Suz */

    misc.ex[i][j] = 0.0;
    if(misc.Suz[i][j] > misc.S[i][j]){
      misc.ex[i][j] = misc.Suz[i][j] - misc.S[i][j];
      misc.Suz[i][j] = misc.S[i][j];
    }

    /* - qv is vertical flow = unsaturated subsurface flow
       - misc.td can be unsaturated zone time delay (td)
         or Effective vertical hydraulic gradient (alpha), depending whether > or < 0
         if td then eq. 18.11; if alpha then eq. 18.12
       - qv limited to the Suz */

    _qv = 0.0;
    if(misc.S[i][j] > 0.0){
      _qv = (params.td > 0.0 ?
	     misc.Suz[i][j] / (misc.S[i][j] * params.td) * params.dt
	     : - params.td * params.K0 * exp(- misc.S[i][j] / params.m));
      if(_qv > misc.Suz[i][j])
	_qv = misc.Suz[i][j];
      misc.Suz[i][j] -= _qv;
      if(misc.Suz[i][j] < ZERO)
	misc.Suz[i][j] = 0.0;
      _qv *= Aatb_r;
    }

    /* sum for total */

    misc.qv[i][j] = _qv;
    misc.qv[i][nidxclass] += misc.qv[i][j];

    /* Calculation of ET -> eq. 18.13
       ET is extracted from the root zone   */

    misc.Ea[i][j] = 0.0;
    if(ET0[i] > 0.0){
      misc.Ea[i][j] = ET0[i] * (1 - misc.Srz[i][j] / params.Srmax);
      if(misc.Ea[i][j] > params.Srmax - misc.Srz[i][j])
	misc.Ea[i][j] = params.Srmax - misc.Srz[i][j];
    }
    misc.Srz[i][j] += misc.Ea[i][j];

    /* Calculation of flow from fully saturated area
       - This section assumes that a/tanB values are ordered from high to low
       - !! mind the difference between the use of idxstats.Aatb_r (1st if) and Aatb_r (2nd if)
       - qo = 0.0 for first topidx class (= upper limit of topidx and fractional area is zero)
       - misc.ex is infiltration excess overland flow 
       - qo is saturation excess overland flow */

    _qo = 0.0;
    if(j > 0){
      if(misc.ex[i][j] > 0.0)
	/* both limits are saturated */
	_qo = idxstats.Aatb_r[j] *
	  (misc.ex[i][j-1] + misc.ex[i][j]) / 2.0;
      else
	if(misc.ex[i][j-1] > 0.0)
	  /* check if lower limit saturated (higher a/tanB value) */
	  _qo = Aatb_r * misc.ex[i][j-1] /
	    (misc.ex[i][j-1] - misc.ex[i][j]) * misc.ex[i][j-1] / 2.0;
    }

    misc.qo[i][j] = _qo;
    misc.qo[i][nidxclass] += misc.qo[i][j];

    misc.ex[i][nidxclass] += Aatb_r * misc.ex[i][j];  /* only for stats */
    misc.qt[i][j] = misc.qo[i][j] + misc.qs[i];	      /* only for stats */
  }


  /* qt = qo + qs */

  /* All dimensions are in [mm] so we can happily add them up.
     The relative area only comes into play at last step.
     fex is independent on the nidxclass, as is qs */

  misc.qo[i][nidxclass] += misc.fex[i];
  misc.qt[i][nidxclass] = misc.qo[i][nidxclass] + misc.qs[i];



  /* S_mean is deficit so qs goes OUT and qv goes in */

  misc.S_mean[i] = misc.S_mean[i] + misc.qs[i] - misc.qv[i][nidxclass];

  if(i + 1 < ntimestep) misc.S_mean[i + 1] = misc.S_mean[i];

  /* qt goes to Qt using the ndelay function */

  for(j=0; j<misc.nreach; j++){
    k = i + j + misc.ndelay;
    if(k > ntimestep - 1)
      break;
    misc.Qt[k] += misc.qt[i][nidxclass] * misc.Ad[j];
  }
	
  return;
}
