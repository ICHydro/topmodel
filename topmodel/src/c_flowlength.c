
/* Calculates the flowlength within the subcatchment area
   Single flow direction */

#include <math.h>
#include <R.h>

void c_flowlength(double *inputdem, double *output, int *nrow, int *ncol, int *iout, int *jout)
{
  int i,j,ii,jj,i2,j2,k,idown,jdown,orig_i,orig_j,flowlength,endflow;
  int *flowi, *flowj;
  double *distance;
  double max,div,divdown,exclude;
  double **flowmap, **dem;
  int no_outlet = 0;

  if(*iout < 0) no_outlet = 1;

  /* memory allocation */

  flowi    = (int *) R_alloc((*nrow) * (*ncol), sizeof(int));
  flowj    = (int *) R_alloc((*nrow) * (*ncol), sizeof(int)); 
  distance = (double *) R_alloc((*nrow) * (*ncol), sizeof(double)); 
                
  dem      = (double **) R_alloc(*nrow, sizeof(double *));
  flowmap  = (double **) R_alloc(*nrow, sizeof(double *));
 
  for(i=0; i<*nrow; i++){
    dem[i]   = (double *) R_alloc(*ncol, sizeof(double));
    flowmap[i]   = (double *) R_alloc(*ncol, sizeof(double));
  }

  /* flowi and flowj store the indexes of the downslope search
     and flowlength counts those */

  exclude = 999999;
  for(k=0; k<*nrow * *ncol; k++){
    flowi[k] = 0;
    flowj[k] = 0;
    distance[k] = 0;
  }

  /* initialise the flowlength map, standard value is -1 (not yet done) */

  for(j=0; j< *ncol; j++){
    for(i=0; i< *nrow; i++){
      dem[i][j] = inputdem[i+(*nrow)*j];
      if(dem[i][j] == exclude) flowmap[i][j] = -2;
      else flowmap[i][j] = -1;
    }
  }

  /* move downstream from each cell */

  orig_i = 0;
  orig_j = 0;
  for(j2=0; j2< *ncol; j2++){
    j = j2;
    for(i2=0; i2< *nrow; i2++){
      i = i2;

      /* If flowmap is -1 we have not yet searched this cell,
         so we start a new flowi and flowj series */

      if(flowmap[i][j] == -1){
        flowlength = 0;
        flowi[0] = i;
        flowj[0] = j;
        distance[0] = 0;
        orig_i = i;
        orig_j = j;
        endflow = 0;
        while(endflow == 0) {
	  /* look for the downstream cell (max downslope dtm difference)
	     If we are on the edge the flowpath ends */
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
                  divdown = div;
                }
              }
            }
          }

	  /* if endflow is 1 here then we've hit a border or an excluded area.
	     If max == 0 we've hit a sink. In all cases the flowpath ends here. */

          if(max == 0) endflow = 1;

          if(endflow == 0) {
            distance[flowlength] = divdown;
            flowlength++;
            flowi[flowlength] = i + idown;
            flowj[flowlength] = j + jdown;
            /* if we arrive in an excluded cell, the flowpath ends too */ 
            if(flowmap[i + idown][j + jdown] == -2) endflow = 1;
          }

          /* if endflow == 1 at this point, we have hit a preliminary end
             and set flowmap to outside (-2), except if no outlet is given */

          if((endflow == 1) && (no_outlet == 0)) {
	    for(k=0; k <= flowlength; k++) {
	      flowmap[flowi[k]][flowj[k]] = -2;
	    }
	  }

	  /* if we reach the outlet, rewrite the distance array inversely
	     to the flowmap */

	  /* Note: iout - 1 and jout -1 because C indices start with 0 */
          
	  /* here is some room for improvement: we don't need to go untill
	     we reach the outlet, once we reach a cell that is alread done we
	     can stop (cfr. river.c) */

	  else {
	    if(((i + idown) == (*iout-1) && (j + jdown) == (*jout-1))
	       || endflow == 1) {
	      if(flowmap[flowi[flowlength]][flowj[flowlength]] < 0) {
		distance[flowlength] = 0;
		flowmap[flowi[flowlength]][flowj[flowlength]] = 0;
	      }
	      else distance[flowlength]
		     = flowmap[flowi[flowlength]][flowj[flowlength]];
	      for(k = (flowlength-1); k >= 0; k--) {
		flowmap[flowi[k]][flowj[k]]
		  = distance[k] + flowmap[flowi[k+1]][flowj[k+1]];
	      }
	      endflow = 1;
	    }
	  }          
	  i = i + idown;
	  j = j + jdown;
	}
	i = orig_i;
	j = orig_j;
      }
    }  
  }

  for(j=0; j < *ncol; j++){
    for(i=0; i < *nrow; i++){
      output[i + (*nrow * j)] = flowmap[i][j];
    }
  }
  return;
}

