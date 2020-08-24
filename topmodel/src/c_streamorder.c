
/* Calculates the streamorder within the subcatchment area
   Single flow direction */

#include <math.h>
#include <R.h>

void c_streamorder(double *inputdem, double *output, int *nrow, int *ncol,
                                                        int *iout, int *jout)
{

  int i,j,ii,jj,i2,j2,idown,jdown,orig_i,orig_j,flowlength;
  int bifurc,order,endflow;
  int flowi[(*nrow) * (*ncol)], flowj[*nrow * *ncol];
  double upslope, downslope, div,exclude, max;
  double **rflag, **dem;

  /* memory allocation */
          
  dem   = (double **) R_alloc(*nrow, sizeof(double *));
  rflag = (double **) R_alloc(*nrow, sizeof(double *));

  for(i=0; i<*nrow; i++)
  {
    dem[i]   = (double *) R_alloc(*ncol, sizeof(double));
    rflag[i]   = (double *) R_alloc(*ncol, sizeof(double));
  }

  /* flowi and flowj store the indexes of the downslope search
     and flowlength counts those */

  exclude = 999999;
  for(i=0; i<*nrow * *ncol; i++){
    flowi[i] = 0;
    flowj[i] = 0;
  }

  /* initialise the rflag map, standard value is -1 (not yet done) */

  for(j=0; j< *ncol; j++){
    for(i=0; i< *nrow; i++){
      dem[i][j] = inputdem[i+(*nrow)*j];
      if(dem[i][j] == exclude) rflag[i][j] = -2;
      else rflag[i][j] = -1;
    }
  }

  /* move downstream from each cell */

  orig_i = 0;
  orig_j = 0;
  for(j2=0; j2< *ncol; j2++){
    j = j2;
    for(i2=0; i2< *nrow; i2++){
      i = i2;

      /* If rflag is -1 we have not yet searched this cell,
         so we start a new flowi and flowj series */

      if(rflag[i][j] == -1){

        /* look for a headwater cell
           there should be no upslope cell, so upslope = 0.0 */

        upslope = 0;
        for(jj=-1; jj < 2; jj++){
          for(ii=-1; ii < 2; ii++){
            if(((i+ii < 0) || (i+ii >= *nrow) || (j+jj < 0) || (j+jj >= *ncol))
               || (dem[i+ii][j+jj] == exclude)) continue;
            if ((ii == 0) && (jj == 0)) continue;
            if((ii == 0) || (jj == 0)) div = 1;
            else div = sqrt(2);
            if(((dem[i+ii][j+jj] - dem[i][j])/div) > max) {
              upslope = ((dem[i+ii][j+jj] - dem[i][j])/div);
            }
          }
        }

          /* but there should be a downslope cell */

        downslope = 0;
        for(jj=-1; jj < 2; jj++){
          for(ii=-1; ii < 2; ii++){
            if(((i+ii < 0) || (i+ii >= *nrow) || (j+jj < 0) || (j+jj >= *ncol))
               || (dem[i+ii][j+jj] == exclude)) continue;
            if ((ii == 0) && (jj == 0)) continue;
            if((ii == 0) || (jj == 0)) div = 1;
            else div = sqrt(2);
            if(((dem[i][j] - dem[i+ii][j+jj])/div) > max) {
              downslope = ((dem[i][j] - dem[i+ii][j+jj])/div);
            }
          }
        }

        /* If this is a headwater cell then we start a flowline */
          
        if((upslope == 0) && (downslope > 0)) {
          flowlength = 1;
          flowi[flowlength-1] = i;
          flowj[flowlength-1] = j;
          order = 1;
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
                  }
                }
              }
            }

            /* if endflow is 1 here we've hit a border or an excluded area.
               If max == 0 we've hit a sink. In all cases the flowpath ends here */

            if(max == 0) endflow = 1;

            if(endflow == 0) {
              flowlength++;
              flowi[flowlength-1] = i + idown;
              flowj[flowlength-1] = j + jdown;
              /* if we arrive in an excluded cell, the flowpath ends too */ 
              if(rflag[i + idown][j + jdown] == -2) endflow = 1;
            }

          /* if endflow = 1 at this point, we have hit a preliminary end
             and set rflag to outside (-2) */

            if(endflow == 1) {
              for(ii=0;ii<flowlength;ii++) {
                rflag[flowi[ii]][flowj[ii]] = -2;
              }
            }

          /* if we arrive at the outlet, trace back and updata the orders
             Note: iout - 1 and jout -1 because C indices start with 0 */

            else {
              if(((i + idown) == (*iout-1)) && (j + jdown) == (*jout-1)) {
                for(ii=0;ii<flowlength;ii++) {
                  /* if we find a cell with a lower order, increase to order */
                  if(rflag[flowi[ii]][flowj[ii]] < order) {
                    rflag[flowi[ii]][flowj[ii]] = order;
                  }
                  else {
                  /* if we find a stream cell from the same order, we need to be
                    careful. It can already be a merger of equals in which case
                    we do not update! */
                    if(rflag[flowi[ii]][flowj[ii]] == order) {
                      bifurc = 0;
                      for(j2=-1; j2 < 2; j2++){
                        for(i2=-1; i2 < 2; i2++){
                          if(((i+i2 < 0) || (i+i2 >= *nrow) 
                            || (j+j2 < 0) || (j+j2 >= *ncol))) continue;
                          if ((i2 == 0) && (j2 == 0)) continue;
                          if(rflag[flowi[ii]+i2][flowj[ii]+j2] == (order-1))
                            bifurc = 1;
                        }
                      }
                    }
                    if(bifurc == 0) {
                      order=order+1;
                      rflag[flowi[ii]][flowj[ii]] = order;
                    }
                  }             
                }
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
      output[i + (*nrow * j)] = rflag[i][j];
    }
  }
  return;
}
