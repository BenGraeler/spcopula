
#include <math.h>
#include <stdio.h>
#include <stdlib.h>



#ifdef _WIN32
#define R_EXPORT __declspec(dllexport)
#else
#define R_EXPORT
#endif




#ifdef __cplusplus
extern "C" {
#endif

/* FORWARD DECLARATIONS */

void R_EXPORT distMatrix(double *out, double *xa, double *xb, int *na, int*nb, int *d);
void R_EXPORT distMatrix_2d(double *out, double *xa, double *xb, int *na, int*nb);
void R_EXPORT distMatrixInner(double *out, double *xa, int *na, int *d);
void R_EXPORT distMatrixInner_2d(double *out, double *xa, int *na);

void R_EXPORT getNeighbours_2d(double *out_locs, double *out_dists, double *data_locs, double *pred_locs, int *n_data, int *n_pred, double *min_dist, int *k, int *prediction);


#ifdef __cplusplus
}
#endif