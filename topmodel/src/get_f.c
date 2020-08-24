#include "topmodel.h"

double
get_f(double t, double R, double C, double K0, double m, double dt)
{

#define	TOLERANCE	0.00001
#define	MAXITER		2000
#define	NTERMS		10

  double        f, f1, f2, fc, R2, sum, diff;
  int           factorial;
  int	        i, j;
  static double cumf, f_, pt, cnst;
  static int    ponding;               /* char OK, could be _Bool or int too */

  /* pt is made a static variable (bug in GRASS?!) to reflect that ponding
     may occur over more than one timestep */

  /* this function can be called more than once, so the static vars
     need to be reset in the beginning of a time series
     (and cannot be initialised during declaration as in GRASS) */
   
  if(t/dt == 1){
    cumf = 0.0;
    f_ = 0.0;
    ponding = 0;
    pt = 0;
  }

  /* R     rainfall intensity (m/h)
     f1    total infiltration (m) at the start of the timestep
     f2    total infiltration (m) at the end of the timestep
     R2    infiltration rate (m/h)
     pt    time to ponding
     cumf  cumulative infiltration (m) at the start of the timestep 
     f_    cumulative infiltration (m) at the end of the timestep
     C    capillary drive, see Morel-Seytoux and Khanji (1974) */

  /* If there is no rainfall go back and reset cumf and ponding */

  if(R <= 0.0){
    cumf = 0.0;
    ponding = 0;
    f_ = 0;
    pt = 0;
    return 0.0;
  }

  f1 = 0.0;
	
  if(!ponding){
    if(cumf > 0){
      f1 = cumf;
      R2 = - K0 / m * (C + f1) / (1 - exp(f1 / m)); /* revise formula!! */
				
      /* if rainfall intensity is higher initial infiltration rate, ponding starts */

      if(R > R2){
	f_ = cumf;
	pt = t - dt;
	ponding = 1;
	goto cont1;
      }
    }

    /* continue if no ponding (yet) */

    f2 = cumf + R * dt;
    R2 = - K0 / m * (C + f2) / (1 - exp(f2 / m));

/* Rprintf("%f, ",R2); */

    /* if rainfall intensity is less than infiltration rate at the end of the timestep.
       All rainfall will infiltrate */

    if(f2 == 0.0 || R < R2){
      f = R;
      cumf += f * dt;
      ponding = 0;
      return f;
    }

    /* Rainfall intensity is greater than infiltration rate at the end of the timestep
       then we need an iterative search for the real cumulative infiltration
       since f2 will be an overestimation and therefore R2 too */

    f_ = cumf + R2 * dt;

    for(i=0; i<MAXITER; i++){
      R2 = - K0 / m * (C + f_) / (1 - exp(f_ / m));
      if(R2 > R){
	f1 = f_;
	f_ = (f_ + f2) / 2.0;
	diff  = f_ - f1;
      }else{
	f2 = f_;
	f_ = (f_ + f1) / 2.0;
	diff  = f_ - f2;
      }
      if(fabs(diff) < TOLERANCE) break;
    }
    /* if the iteration failed return -9999 */
		
    if(i == MAXITER){
      f = -9999;
      return f;
    }

    /* the iteration may result in no ponding after all */
		
    pt = t - dt + (f_ - cumf) / R;
    if(pt > t){
      f = R;
      cumf += f * dt;
      ponding = 0;
      return f;
    }

    /* if we make it up to here there is ponding.
       At this point f_ is cumulative infiltration
       at start of ponding (Ip in Beven 1984)
       so we need to add infiltration during ponding */

    /* calculation of cnst (lambda in Beven (1984)) is done
       when ponding is reached and is kept until ponding stops */

  cont1:
    cnst = 0.0;
    factorial = 1;
    fc = (f_ + C);
    for(j=1; j<=NTERMS; j++){
      factorial *= j;
      cnst += pow(fc / m, (double) j) /
	(double) (j * factorial);
    }
    cnst = log(fc) - (log(fc) + cnst) / exp(C / m);
    f_ += R * (t - pt) / 2.0;         /* first estimation of f_  */
    ponding = 1;

  }

  /* if ponding occurs then we need to calculate the ponding infiltration rate.
     (eq. (6) in beven 1984). This is easy with a Newton-Raphson iteration */

  /* fc = f_ + C is used to make calculation easier to read */
  /* cnst = lambda in Beven, 1984 */

  /* Newton Raphson: f1 = f(x); f2 = f'(x) and f = - f(x)/f'(x)
     (see e.g. wikipedia) */  

  for(i=0; i<MAXITER; i++){
    fc  = f_ + C;
    sum = 0.0;
    factorial = 1;
    for(j=1; j<=NTERMS; j++){
      factorial *= j;
      sum += pow(fc / m, (double) j) / (double) (j * factorial);
    }
    f1  = - (log(fc) - (log(fc) + sum) / exp(C / m) - cnst) /
      (K0 / m) - (t - pt);
    f2  = (exp(f_ / m) - 1.0) / (fc * K0 / m);
    diff   = - f1 / f2;
    f_ += diff;
    if(fabs(diff) < TOLERANCE) break;
  }

  if(i == MAXITER){
    f = -9999;
    return f;
  }

  if(f_ - cumf < R * dt){	               /* only in this case ponding */
    f    = (f_ - cumf) / dt;
    cumf = f_;
    /* initial guess for next timestep */
    f_  += f * dt;
    return f;
  }else{
    f = R;
    cumf += f * dt;
    ponding = 0;
    pt = 0;
    return f;
  }
}

