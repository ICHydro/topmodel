
/* filling sinks in a dem */

#include <math.h>
#include <R.h>

void c_sinkfill(double *input, double *output, int *nrow, int *ncol, double *cellsize, double *degree)
{

  int	numsin, i, j, xloc, yloc, sinkflag;
  int inix, iniy,itermax, numiter, numoutside;
  double grad, elemin, inimin, outside, degree2;
  double **dem;

  /* memory allocation */
  	
  dem = (double **) R_alloc(*nrow, sizeof(double *));
  for(i=0; i<*nrow; i++){
    dem[i] = (double *) R_alloc(*ncol, sizeof(double));
  }

  /* copy input to dem (R fills matrices per column) */

  for(j=0; j< *ncol; j++){
    for(i=0; i< *nrow; i++){
      dem[i][j] = input[i+(*nrow)*j];
    }
  }

  /* initialisation of other variables */

  outside = 999999;
  itermax = 100;
  degree2 = *degree / 360 * (2 * 3.1415927);
  inimin = 2000000;
  numiter = 0;
  numsin = -9999;
  inix = 0;
  iniy = 0;

  /* find the outlet = cell with minimum elevation */

  for(j = 0; j < *ncol; j++) {
    for(i = 0; i < *nrow; i++) {
      if(dem[i][j] < -1000) dem[i][j] = 999999;
      if(dem[i][j] < inimin) {
	inimin = dem[i][j];
	inix = i;
	iniy = j;
      }
    }
  }

  /*  main loop for sink filling */

  while((numsin != 0) && (numiter < itermax)) {

    numiter = numiter + 1;
    numsin = 0;

    for(j = 1; j < (*ncol-1); j++) {
      for(i = 1; i < (*nrow-1); i++) {
	sinkflag = 1;
	numoutside = 0;
	if((dem[i][j] <= 0) || (dem[i][j] == outside) || ((i == inix) && (j == iniy)))
	  continue;               /* skip cells outside and the outlet */
	for(yloc = -1; yloc < 2; yloc++) {
	  for(xloc = -1; xloc < 2; xloc++) {
	    if((xloc != 0) || (yloc != 0)) {
	      if((xloc == 0) || (yloc == 0)) grad = tan(degree2) * *cellsize;
	      else grad = tan(degree2) * sqrt(pow(*cellsize,2) + pow(*cellsize,2));
	      if(dem[i][j] > dem[i+xloc][j+yloc] + grad - 0.00000001) sinkflag = 0;
	      if((dem[i+xloc][j+yloc] == outside) && numiter > 10) numoutside++;	/* why numiter > 10 ???? */	
	    }
	  }
	}

	/*  if sinkflag ==1, then we have a sink except if it borders an outside cell */

	if((sinkflag == 1) && (numoutside == 0)) {
	  numsin++;
	  elemin = 999999999;

	  /*  look for the minimum elevation around a cell and increase the elevation of sink cell */

	  for(yloc = -1; yloc < 2; yloc++) {
	    for(xloc = -1; xloc < 2; xloc++) {
	      if((dem[i+xloc][j+yloc] != outside) && ((xloc != 0) || (yloc != 0))) {
		if(dem[i+xloc][j+yloc] < elemin) {
		  elemin = dem[i+xloc][j+yloc];
		  if((xloc == 0) || (yloc == 0)) grad = tan(degree2) * *cellsize;
		  else grad = tan(degree2) * sqrt(pow(*cellsize,2) + pow(*cellsize,2));
		}
	      }
	    }
	  }
	  dem[i][j] = elemin + grad;
	}
      }
    }
  }

  /*  format the output */

  output[0] = numiter;
  output[1] = numsin;

  /* R fills matrices per column */

  for(j=0; j < *ncol; j++){
    for(i=0; i < *nrow; i++){
      output[2 + i + (*nrow * j)] = dem[i][j];
    }
  }
	
  return;
}

