#ifndef SPARSE_MATRIX_H
#define SPARSE_MATRIX_H

#include <algorithm>
#include <gsl/gsl_spmatrix.h>

struct SparseMat_struct {
  int nrow,ncol;

  gsl_spmatrix *gslmatrix;

  int nz, nzmax;
  int *indx, *indy;
  double *data, tol;

  /* constructor */
  SparseMat_struct(int N, int P);

  /* destructor */
  ~SparseMat_struct();

  /* methods */
  double get(int i, int j);
  void set(int i, int j, double val);
  void set_zero();

  void init();
  void extend_storage(int newnzmax);
};

typedef SparseMat_struct SparseMat;

#endif
