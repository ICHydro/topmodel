/* Find rivers based on some thresholds of topidx and area

   Defining an appropriate mask automatically defines the outlet
   (lowest point at the border of the catchment).

   If no mask present, the rivers stop at the border. The algorithm finds
   all headwater cells, traces down, and remaps the distance from the outlet.

   The "river" matrix is coded as follows:
   -999999: outside catchment
   -10: no river (default value, initialised)
   -1:  first initialisation as river cell
   -2:  identified as headwater cell
   >=0:  river cell with calculated distance */

#include <math.h>
#include <R.h>

void findrivers(double *inputdem, double *inputtopidx, double *inputarea,
		double *outputriver, int *nrow, int *ncol, double *cellsize,
		double *thtopidx, double *tharea)
{

  int i,j,k,idown,jdown,river_origin_i,river_origin_j,endriver,nriv;
  int ii,jj,i2,j2,headwater;
  int *riveri, *riverj;

  double max, addist, exclude, dist3;
  double **dem, **atb, **area, **river, *dist;

  /* memory allocation */

  riveri= (int *) R_alloc((*nrow) * (*ncol), sizeof(int));
  riverj= (int *) R_alloc((*nrow) * (*ncol), sizeof(int));
  dist  = (double *) R_alloc((*nrow) * (*ncol), sizeof(double));
          
  dem   = (double **) R_alloc(*nrow, sizeof(double *));
  atb   = (double **) R_alloc(*nrow, sizeof(double *));
  area  = (double **) R_alloc(*nrow, sizeof(double *));
  river = (double **) R_alloc(*nrow, sizeof(double *));

  for(i = 0; i < *nrow; i++){
    dem[i]   = (double *) R_alloc(*ncol, sizeof(double));
    atb[i]   = (double *) R_alloc(*ncol, sizeof(double));
    area[i]  = (double *) R_alloc(*ncol, sizeof(double));
    river[i] = (double *) R_alloc(*ncol, sizeof(double));
  }

  /* copy input to dem (R fills matrices per column) */

  for(j=0; j< *ncol; j++){
    for(i=0; i< *nrow; i++){
      dem[i][j] = inputdem[i+(*nrow)*j];
      atb[i][j] = inputtopidx[i+(*nrow)*j];
      area[i][j] = inputarea[i+(*nrow)*j];
    }
  }

  /* initialisation of other variables */

  exclude = -9999;
 
  /* identify river cells */

  for(j=0; j< *ncol; j++){
    for(i=0; i< *nrow; i++){
      river[i][j] = -10;
      if(dem[i][j] == exclude) river[i][j] = -999999;
      if((area[i][j] > *tharea) && (atb[i][j] > *thtopidx) && dem[i][j] > exclude)
	river[i][j] = -1;
    }
  }

  /* search for headwater cells (all cells with no other river cells upslope) */
 
  for(j=0; j< *ncol; j++){
    for(i=0; i< *nrow; i++){
      headwater = 1;
      if(river[i][j] != -1) continue;
      for(jj=-1; jj < 2; jj++){
	for(ii=-1; ii < 2; ii++){
          if(((i+ii >= 0) && (i+ii < *nrow) && (j+jj >= 0) && (j+jj < *ncol))
             && ((ii != 0) || (jj != 0)) 
             && (dem[i+ii][j+jj] != exclude)
             && ((river[i+ii][j+jj] == -1) || (river[i+ii][j+jj] == -2))
             && (dem[i][j] < dem[i+ii][j+jj]))
            headwater = 0;           
        }
      }
      if(headwater == 1) river[i][j] = -2;
    }
  }

  /* follow the rivers downstream
     i2 and j2 are used for the main loop,
     i and j follow the river downstream
     and are reset to i2 and j2 at the end of the river tracking */

  for(j2=0; j2< *ncol; j2++){
    j = j2;
    for(i2=0; i2< *nrow; i2++){
      i = i2;
      if((dem[i][j] == exclude) || river[i][j] != -2) continue;
      river[i][j] = -1;
      endriver = 0;
      nriv = 1;
      riveri[nriv] = i;
      riverj[nriv] = j;
      dist[1] = 0.0;
      river_origin_i = i;
      river_origin_j = j;
      while(endriver == 0) {
        max = 0.0;
        for(jj=-1; jj < 2; jj++){
          for(ii=-1; ii < 2; ii++){
            if(((i+ii >= 0) && (i+ii < *nrow) && (j+jj >= 0) && (j+jj < *ncol))
               && ((ii != 0) || (jj != 0)) 
               && (dem[i+ii][j+jj] != exclude)) {
              if((ii == 0) || (jj == 0)) dist3 = *cellsize;
              else dist3 = sqrt(pow(*cellsize,2) + pow(*cellsize,2));
              if(((dem[i][j] - dem[i+ii][j+jj])/dist3) > max) {
                max = ((dem[i][j] - dem[i+ii][j+jj])/dist3);
                idown = ii;
                jdown = jj;
                addist = dist3;
              }
            }
          }
        }

        /* idown and jdown now give the relative position of the dowslope cell
           addist gives the distance, max gives the difference in altitude

           if max = 0, then we are at a sink or at an outlet
           (either at the border of the grid or the exclude area): */

        if(max == 0) endriver = 1;

        if(endriver != 1) {
          dist[nriv] = addist;
          nriv = nriv + 1;
          riveri[nriv] = i + idown;
          riverj[nriv] = j + jdown;

          /* if we pass a headwater cell, reset this to a normal river cell */

	  if(river[i+idown][j+jdown] == 2) river[i+idown][j+jdown] = -1;

          /*  two more conditions for finishing:
              - trace flows into an existing river
              - river ends outside the catchment */

          if(river[i+idown][j+jdown] == 0) endriver = 1;
          if(river[i+idown][j+jdown] == exclude) endriver = 1;
        }

        /* if at the end of the river, rewrite the distance array inversely to rivermap */

        if(endriver == 1) {
          if(river[riveri[nriv]][riverj[nriv]] < 0) {
            dist[nriv] = 0;
            river[riveri[nriv]][riverj[nriv]] = 0;
          }
          else dist[nriv] = river[riveri[nriv]][riverj[nriv]];
          for(k = (nriv - 1); k > 0; k--) {
            river[riveri[k]][riverj[k]] = dist[k] + river[riveri[k+1]][riverj[k+1]];
          }
        }
        else {
          i = i + idown;
          j = j + jdown;
        }
      }
      i = river_origin_i;
      j = river_origin_j;
    }
  }

  for(j=0; j < *ncol; j++){
    for(i=0; i < *nrow; i++){
      outputriver[i + (*nrow * j)] = river[i][j];
    }
  }
  return;
}
