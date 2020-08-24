
/* Delineates the subcatchment area based on outlet coordinates and
   a single flow direction algorithm */

#include <math.h>
#include <R.h>

void c_subcatch(double *inputdem, int *output, int *nrow, int *ncol, int *iout, int *jout)
{
  int i,j,ii,jj,exclude,i2,j2,endflow;
  int idown,jdown,origin_i,origin_j,flowlength;
  int *flowi, *flowj;
  int **catchment;
  double **dem;
  double max,div;

  /* flowi and flowj store the indexes of the downslope flow search */

  /* memory allocation */
  	
  dem        = (double **) R_alloc(*nrow, sizeof(double *));
  catchment  = (int **) R_alloc(*nrow, sizeof(int *));
  flowi      = (int *) R_alloc((*nrow) * (*ncol), sizeof(int));
  flowj      = (int *) R_alloc((*nrow) * (*ncol), sizeof(int));

  for(i=0; i<*nrow; i++){
    dem[i]        = (double *) R_alloc(*ncol, sizeof(double));
    catchment[i]  = (int *) R_alloc(*ncol, sizeof(int));
  }
  
  /* copy input to dem (R fills matrices per column) */

  for(j=0; j< *ncol; j++){
    for(i=0; i< *nrow; i++){
      dem[i][j] = inputdem[i+(*nrow)*j];
    }
  }

  /* initialisation */

  exclude = 999999;
  flowlength = 0;
  max = 0;

  for(i=0; i < *ncol * *nrow; i++) {
    flowi[i] = 0;
    flowj[i] = 0;
  }

  /* initialise the catchment map, standard value is -1
     (outside the catchment value is 0) */

  for(j=0; j< *ncol; j++){
    for(i=0; i< *nrow; i++){
      catchment[i][j] = -1;
      if(dem[i][j] == exclude) catchment[i][j] = 0;
    }
  }

  /* NOTE: iout-1 and jout-1 because R indices start at 1 while at 0 in C */

  catchment[*iout-1][*jout-1] = 1;

  /* move downstream from each cell */

  for(j2=0; j2< *ncol; j2++){
    j = j2;
    for(i2=0; i2< *nrow; i2++){
      i = i2;

      /* If catchment is -1 we have not yet searched this cell,
	 so we start a new flowi and flowj series */

      if(catchment[i][j] == -1) {
	flowlength = 1;
	flowi[flowlength-1] = i;
	flowj[flowlength-1] = j;
	origin_i = i;
	origin_j = j;
	endflow = 0;
	while(endflow == 0) {

	  /* look for the downstream cell (max downslope dem difference) */

	  max = 0;
	  for(jj=-1; jj < 2; jj++){
	    for(ii=-1; ii < 2; ii++){
	      if(((i+ii < 0) || (i+ii >= *nrow) || (j+jj < 0) || (j+jj >= *ncol))
		 || (dem[i+ii][j+jj] == exclude)) endflow = 1;
	      else {
		if ((ii == 0) && (jj == 0)) continue;
		if((ii == 0) || (jj == 0)) div = 1;
		else div = sqrt(2);
		if(((dem[i][j] - dem[i+ii][j+jj])/div) > max) {
 		  max = ((dem[i][j] - dem[i+ii][j+jj])/div);
		  idown = ii;
		  jdown = jj;
		}
	      }
	    }
	  }

	  /* if endriver is 1 at this moment we've hit a border or an excluded area
	     in both cases the flowpath is not included in the catchment */

	  /* if max ==0 we've hit a sink, in which case the trace is also disregarded */

          if(max == 0) endflow = 1;

	  if(endflow != 1) {
	    flowlength++;
	    flowi[flowlength-1] = i + idown;
	    flowj[flowlength-1] = j + jdown;
            if(catchment[i+idown][j+jdown] == 0) endflow = 1;
	  }

	  if(endflow == 1) {
	    for(ii=0; ii < flowlength; ii++) {
	      catchment[flowi[ii]][flowj[ii]] = 0;
	    }
	  }
	  else {

	    /* if we arrive at the outlet, set all the catchments
	       of the trace to inside (1) */

	    if(catchment[i+idown][j+jdown] == 1) {
	      for(ii=0; ii < flowlength; ii++) {
		catchment[flowi[ii]][flowj[ii]] = 1;
	      }
	      endflow = 1;
	    }
	  }
	  i = i + idown;
	  j = j + jdown;
	}
	i = origin_i;
	j = origin_j;
      }
    }
  }

  for(j=0; j < *ncol; j++){
    for(i=0; i < *nrow; i++){
      output[i + (*nrow * j)] = catchment[i][j];
    }
  }
  return;
}
