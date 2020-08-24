
/* calculates the topographic index (a/ tan(beta) ) */

/*  Average downslope slope is also calculated but not returned */

#include <math.h>
#include <R.h>

#define	ZERO		0.0000001

void c_topidx(double *inputdem,
	    int    *inputriver,
	    int    *nrow, 
	    int    *ncol,
	    double *ew_res,
	    double *ns_res,
	    double *output)
{
  int    i,j,ii,jj,k,k1,k2,k3,nf,i2,j2,natbold;
  int    natb = 0;     /* number of atbs still to be analysed */
  int    iproute,nrout,river,not_yet;
  double **dem, **atb, **area, **slope, **rivermap;
  double exclude,dnx,routefac,nslp,xr,yr;
  int	 nsink = 0;
  double routdem[9], tanb[9];
  double c,dx1,dx2,sum,sumtb;

  /* memory allocation */
          
  dem   = (double **) R_alloc(*nrow, sizeof(double *));
  atb   = (double **) R_alloc(*nrow, sizeof(double *));
  area  = (double **) R_alloc(*nrow, sizeof(double *));
  slope = (double **) R_alloc(*nrow, sizeof(double *));
  rivermap = (double **) R_alloc(*nrow, sizeof(double *));

  for(i=0; i<*nrow; i++){
    dem[i]   = (double *) R_alloc(*ncol, sizeof(double));
    atb[i]   = (double *) R_alloc(*ncol, sizeof(double));
    area[i]  = (double *) R_alloc(*ncol, sizeof(double));
    slope[i] = (double *) R_alloc(*ncol, sizeof(double));
    rivermap[i] = (double *) R_alloc(*ncol, sizeof(double));
  }

  /* copy input to dem (R fills matrices per column) */

  for(j=0; j< *ncol; j++){
    for(i=0; i< *nrow; i++){
      dem[i][j] = inputdem[i+(*nrow)*j];
      rivermap[i][j] = inputriver[i+(*nrow)*j];
    }
  }

  /* initialisation */
  /* atb initalised with elevation, -9999 for excluded cells */

  exclude = -9999;
  river = 0;
  dx1 = 1 / *ew_res;
  dx2 = 1 / (sqrt(2) * *ew_res);
  natbold = -1;

  for(j=0; j< *ncol; j++){
    for(i=0; i< *nrow; i++){
      if(dem[i][j] < -9000) dem[i][j] = exclude;
      area[i][j] = *ew_res * *ns_res;
      if(dem[i][j] == exclude) {
        natb++;
        atb[i][j] = exclude;
      } else {
        atb[i][j] = -9.9;
      }
    }
  }

  /* if natbold = natb after the loop then the loop did not fill any more cells
     and we are probably in an indefinite loop so quit.
     The unfinished areas will show up as NA */

  /* an alternative is to work with a fixed loop, and some output to the user:

     Rprintf("ncells = %i\n", ncells);
     Rprintf("natb = %i\n", natb);

     Rprintf("Iterations:         ");

     for(iter=1;natb<ncells;iter++){

     R_CheckUserInterrupt();
     Rprintf("\b\b\b\b\b\b\b\b%8i",iter); */

  while((natb != natbold) && (natb < *nrow * *ncol)) {

    natbold = natb;

    /* loop through the grid cell and check if there is an upslope element
       without an atanb value. If so, we cannot yet calculate the topidx. */

    for(j = 0; j < *ncol; j++) {
      for(i = 0; i < *nrow; i++) {
              
        /* skip non catchment cells and cells that are done */
        if((dem[i][j] == exclude) || (atb[i][j] >= ZERO))
	  continue;

        /* river cells don't accumulate flow downstream and use the average of
           the inflow slope as the local gradient */

        if(rivermap[i][j] == 1) river = 1;
        else {

          /* check the 8 flow directions for upslope elements 
             without a topidx value */

          not_yet = 0;

          for(jj=-1; jj < 2; jj++){
            for(ii=-1; ii < 2; ii++){
              if(((i+ii >= 0) && (i+ii < *nrow) && (j+jj >= 0) && (j+jj < *ncol))
                 && ((ii != 0) || (jj != 0)) 
                 && (dem[i+ii][j+jj] != exclude)) {
                if((dem[i+ii][j+jj] > dem[i][j]) && (atb[i+ii][j+jj] < ZERO))
                  not_yet = 1;
              }
            }
          }

          if(not_yet) continue;

          /* if there are no upslope elements without a topidx value,
             start calculations */

          /* find the outflow direction and calculate the sum of weights using 
             (tanb*countour length). Contour length = 0.5dx for the cardinal
             direction and 0.354dx for diagonal */

          sum = 0;
          for(ii=0; ii<9; ii++){
            tanb[ii] = 0;
            routdem[ii] = 0;
          }

          k = 0;
          sumtb = 0;
          nrout = 0;
          for(jj=-1; jj < 2; jj++){
            for(ii=-1; ii < 2; ii++){
              if(((i+ii >= 0) && (i+ii < *nrow) && (j+jj >= 0) && (j+jj < *ncol))
                 && ((ii != 0) || (jj != 0)) 
                 && (dem[i+ii][j+jj] != exclude)) {
                if((ii == 0) || (jj == 0)) {
                  dnx = dx1;
                  routefac = 0.5;
                } else {
                  dnx = dx2;
                  routefac = 0.5 / sqrt(2);
                }
                if(dem[i][j] - dem[i+ii][j+jj] > ZERO ) {
                  tanb[k] = (dem[i][j] - dem[i+ii][j+jj]) * dnx;
                  routdem[k] = routefac * *ew_res * tanb[k];
                  sum += routdem[k];
                  sumtb = sumtb + tanb[k];
                  nrout++;
                }
              }
              k++;
            }
          }
        }

	/* label 1 */

	/* if a sink or a river cell... */

        if((nrout == 0) || (river == 1)) {
	  nsink++;
          river = 0;

          /* assume that there is a channel of length dx running midway through
             the sink or boundary node. Take average inflow slope angle to
             represent tanb and A/(2dx) to represent a */

          sumtb = 0;
          nslp = 0;
          for(jj=-1; jj < 2; jj++){
            for(ii=-1; ii < 2; ii++){
              if(((i+ii >= 0) && (i+ii < *nrow) && (j+jj >= 0) && (j+jj < *ncol))
                 && ((ii != 0) || (jj != 0)) 
                 && (dem[i+ii][j+jj] != exclude)) {
                if((ii == 0) || (jj == 0)) dnx = dx1;
                else dnx = dx2;
                if(rivermap[i][j] == 0) {
                  /* for sink/boundary squares sumtb is just the average slope */
                  sumtb += (dem[i+ii][j+jj] - dem[i][j]) * dnx;
                  nslp++;
                } else {
                  /* for a river cell the slope is the average inflow slope */
                  if(dem[i+ii][j+jj] > dem[i][j]) {
                    sumtb += (dem[i+ii][j+jj] - dem[i][j]) * dnx;
                    nslp++;
                  }
                }
              }
            }
          }

          /* calculate the average inflow slope angle */

          if(sumtb > ZERO) {
            sumtb = sumtb / nslp;
            atb[i][j] = log(area[i][j] / (2 * sumtb));
            slope[i][j] = 2 * sumtb;
          } else {
            atb[i][j] = exclude;
          }
          natb++;
          continue;
        }

        /*  so much for rivers and sinks. Go on with normal river cells
            NOTE: no mechanism in place to remove negative atb values
            Can be done in postprocessing */

        if(sum > 0) {
          c = area[i][j] / sum;
          atb[i][j] = log(c);
          slope[i][j] = sumtb / nrout;
        } else {
          c = 0;
          atb[i][j] = 0.01;
          slope[i][j] = 0;
        }
        natb++;

        /* calculate downslope area */

	nrout = 0;
	for(jj=-1; jj < 2; jj++){
	  for(ii=-1; ii < 2; ii++){
	    if(((i+ii >= 0) && (i+ii < *nrow) && (j+jj >= 0) && (j+jj < *ncol))
	       && ((ii != 0) || (jj != 0)) 
	       && (atb[i+ii][j+jj] != exclude)) {
	      if(routdem[nrout] > 0) area[i+ii][j+jj] += c * routdem[nrout];
	    }
	    nrout++;
	  }
	}
      }
    }
  }

  /*  Rprintf("\nNumber of sinks or boundaries: %i\n",nsink); */

  /* format output */

  for(j=0; j < *ncol; j++){
    for(i=0; i < *nrow; i++){
      if(atb[i][j] == exclude) area[i][j] = exclude;
      output[i + (*nrow * j)] = atb[i][j];
      output[*nrow * *ncol + i + (*nrow * j)] = area[i][j];
    }
  }
  return;
}




