
#include "Rinterface.h"





/* Helper Functions */
inline double dist_2d(double ax, double ay, double bx, double by) {
	return(sqrt( (ax-bx)*(ax-bx) + (ay-by)*(ay-by) ));
}



#ifdef __cplusplus
extern "C" {
#endif






/**
* Computes the distance matrix of two sets of n-dimensional points
* @param out result (na,nb) matrix (mast be allocated in R or in calling function before)
* @param xa coordinates of the 1st set of points, each row of this matrix is a point
* @param xb coordinates of the 2nd set of points, each row of this matrix is a point
* @param na number of points in xa
* @param nb number of points in xb
* @param d  number of dimensions, equals the number of columns in xa and xb
**/ 
void R_EXPORT distMatrix(double *out, double *xa, double *xb, int *na, int*nb, int *d) {
	if (*d == 2) {
		distMatrix_2d(out, xa, xb, na, nb);
	}
	else
	{
		for (int i=0;i<*na; ++i) {
			for (int j=0;j<*nb; ++j) {
				out[j * (*na) + i] = 0;
				
				for (int k=0; k<*d; ++k) {
					out[j * (*na) + i] += (xa[k * (*na) + i] -  xb[k * (*nb) + j]) * (xa[k * (*na) + i] -  xb[k * (*nb) + j]); // In R: col major
					
				}
				out[j * (*na) + i] = sqrt(out[j * (*na) + i]);		
			}
		}
	}
}



/**
* Computes the distance matrix of two sets of 2-dimensional points
* @param out result (na,nb) matrix (mast be allocated in R or in calling function before)
* @param xa coordinates of the 1st set of points, each row of this matrix is a point
* @param xb coordinates of the 2nd set of points, each row of this matrix is a point
* @param na number of points in xa
* @param nb number of points in xb
**/ 
void R_EXPORT distMatrix_2d(double *out, double *xa, double *xb, int *na, int*nb) {
	for (int i=0;i<*na; ++i) {
		for (int j=0;j<*nb; ++j) {
			out[j * (*na) + i] =  dist_2d(xa[i],xa[*na + i], xb[j], xb[*nb + j]);	
		}
	}
}



/**
* Computes the "inner" distance matrix of one set of n-dimensional points
* @param out result (na,na) matrix (mast be allocated in R or in calling function before)
* @param xa coordinates of the 1st set of points, each row of this matrix is a point
* @param na number of points in xa
* @param d  number of dimensions, equals the number of columns in xa and xb
**/ 
void R_EXPORT distMatrixInner(double *out, double *xa, int *na, int *d) {
	if (*d == 2) {
		distMatrixInner_2d(out, xa, na);
	}
	else
	{
		for (int i=0;i<*na; ++i) {
			for (int j=i+1;j<*na; ++j) {
				out[j * (*na) + i] = 0;
				for (int k=0; k<*d; ++k) {
					out[j * (*na) + i] += (xa[k * (*na) + i] -  xa[k * (*na) + j]) * (xa[k * (*na) + i] -  xa[k * (*na) + j]); // In R: col major
					
				}
				out[j * (*na) + i] = sqrt(out[j * (*na) + i]);	
				out[i * (*na) + j] = out[j * (*na) + i]; 
			}
		}
		for (int i=0;i<*na; ++i) {
			out[i * (*na) + i] = 0.0;
		}
	}
}



/**
* Computes the "inner" distance matrix of one set of n-dimensional points
* @param out result (na,na) matrix (mast be allocated in R or in calling function before)
* @param xa coordinates of the 1st set of points, each row of this matrix is a point
* @param na number of points in xa
**/ 
void R_EXPORT distMatrixInner_2d(double *out, double *xa, int *na) {
	for (int i=0;i<*na; ++i) {
		for (int j=i+1;j<*na; ++j) {
			out[j * (*na) + i] =  dist_2d(xa[i],xa[*na + i], xa[j], xa[*na + j]);	
			out[i * (*na) + j] = out[j * (*na) + i]; 
		}
	}
	// set diagonal to zero
	for (int i=0;i<*na; ++i) {
		out[i * (*na) + i] = 0.0;
	}
}






/*for (i in 1:nLocs) {
    tempDists <- spDists(dataLocs, predLocs[i, ])
    tempDists[tempDists < min.dist] <- Inf
    spLocs <- order(tempDists)[1:(size - 1)]
    
    allLocs[i,] <- c(i, spLocs)
    allDists[i,] <- tempDists[spLocs]
    allData[i,(prediction+1):size] <- dataLocs[c(i[!prediction],spLocs),
                                               var, drop = F]@data[[1]]
  }
  */

/**
*
*
*
*
* @param k size of the neighbourhood including the location of interest
*/
void R_EXPORT getNeighbours_2d(double *out_locs, double *out_dists, double *data_locs, double *pred_locs, int *n_data, int *n_pred, double *min_dist, int *k, int *prediction) {
	double *C = (double*)malloc((*n_data) *  (*n_pred) * sizeof(double));
	double *knn_dists = (double*)malloc((*k-1) * sizeof(double));
	int *knn_idxs  = (int*)malloc((*k-1) * sizeof(int));
	
	// Depending on whether pred_locs are the same as data_locs, use optimized "inner" distance matrix calculation
	if (!prediction) { 
		// Compute inner distance matrix and store in C
		distMatrixInner_2d(C, data_locs, n_data);
	}
	else {
	
		// Compute distance matrix and store in C
		distMatrix_2d(C, data_locs, pred_locs, n_data, n_pred);
	}	
	
	// for each pred location
	for (int i=0; i<*n_pred; ++i) {
		knn_dists[0] = INFINITY;
		knn_idxs[0] = -1;
		
		for (int j=0; j<*n_data; ++j) {
			if (C[i* (*n_data) + j] < *min_dist) continue;
			int nbrank = 0;
			while (C[i* (*n_data) + j] > knn_dists[nbrank] &&  nbrank < (*k-1)) ++nbrank; 
			if (nbrank < (*k-1)) {		
				for (int s=(*k-3); s >= nbrank; --s) {
					knn_dists[s + 1] = knn_dists[s]; // shift array if new element in knn array
					knn_idxs[s + 1] = knn_idxs[s];
				}
				knn_dists[nbrank] = C[i* (*n_data) + j];
				knn_idxs[nbrank] = j+1; // Indices in R start with 1
			}			
		}
		
		out_locs[i] = i+1; // Indices in R start with 1
		for (int j=1; j<(*k); ++j) {			
			// Assign c(i, knn_idxs[1:size-1]) as i-th row in out_locs
			out_locs[j*(*n_pred) + i] = knn_idxs[j-1];			
			// Assign knn_dists[1:size-1] as i-th row in out_dists
			out_dists[(j-1)*(*n_pred) + i] = knn_dists[j-1];
		}
		
	}
	
	free(knn_idxs);
	free(knn_dists);
	free(C);
}






/*
calcSpLagInd <- function(data, boundaries) {
  lags <- vector("list",length(boundaries))
  
  dists <- spDists(data)
  nlocs <- length(data)
  
  for (i in 1:(nlocs-1)) {
    for (j in (i+1):nlocs) {
      d <- dists[i,j]
      for ( k in 1:length(boundaries)) {
        if (d < boundaries[k]) {
          lags[[k]] <- rbind(lags[[k]],c(i,j,d))
          break()
        }
      }
    }
  }
  return(lags)
}
*/

/**
*
* @
*
*
*
* OpenMP?
*/
// SEXP calcSpLagInd(double *o double *idata, double *iboundaries, int *in_boundaries) {


	// for

// }



#ifdef __cplusplus
}
#endif