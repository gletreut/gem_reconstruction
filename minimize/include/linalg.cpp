//****************************************************************************
// linalg.cpp - Implementation of the definitions in linalg.h.
// Date: 2017-05-03
// Created by: Guillaume Le Treut
//****************************************************************************

//#define GSL
//#define LAPACK
//#define MKL

/* custom files */
#include "linalg.h"

using namespace std;
/* linear algebra */
double
linalg_dnrm2(const gsl_vector *v){
  double vnorm;
#ifdef GSL
  vnorm = gsl_blas_dnrm2(v);
#elif defined(LAPACK) || defined(MKL)
  vnorm = cblas_dnrm2(v->size, v->data, 1);
#endif

  return vnorm;
}

double
linalg_dnrm2(const gsl_matrix *m){
  gsl_vector_view mview = gsl_vector_view_array(m->block->data,m->block->size);

  return linalg_dnrm2(&mview.vector);
}

double
linalg_ddot(const gsl_vector *a, const gsl_vector *b){
  double ddot(0);

#ifdef GSL
  gsl_blas_ddot(a,b,&ddot);
#elif defined(LAPACK) || defined(MKL)
  ddot = cblas_ddot(a->size, a->data, 1, b->data, 1);
#endif
  return ddot;
}

double
linalg_ddot(const gsl_matrix *a, const gsl_matrix *b){
  gsl_vector_view aview = gsl_vector_view_array(a->block->data,a->block->size);
  gsl_vector_view bview = gsl_vector_view_array(b->block->data,b->block->size);

  return linalg_ddot(&aview.vector, &bview.vector);
}

void
linalg_dscal(const double s, gsl_vector *a){
#ifdef GSL
  gsl_blas_dscal(s,a);
#elif defined(LAPACK) || defined(MKL)
  cblas_dscal(a->size, s, a->data, 1);
#endif
  return;
}

void
linalg_dscal(const double s, gsl_matrix *a){

  gsl_vector_view aview = gsl_vector_view_array(a->block->data,a->block->size);

  linalg_dscal(s,&aview.vector);
  return;
}

void
linalg_daxpy(const double alpha, const gsl_vector *a, gsl_vector *b){

#ifdef GSL
  gsl_blas_daxpy(alpha,a,b);
#elif defined(LAPACK) || defined(MKL)
  cblas_daxpy(a->size, alpha, a->data, 1, b->data, 1);
#endif
  return;
}

void
linalg_daxpy(const double alpha, const gsl_matrix *a, gsl_matrix *b){

  gsl_vector_view aview = gsl_vector_view_array(a->block->data,a->block->size);
  gsl_vector_view bview = gsl_vector_view_array(b->block->data,b->block->size);

  linalg_daxpy(alpha,&aview.vector,&bview.vector);
  return;
}

void
linalg_dgemm(int transA, int transB, const double alpha, const gsl_matrix *A, const gsl_matrix *B, const double beta, gsl_matrix *C){

#ifdef GSL
  CBLAS_TRANSPOSE_t transA_, transB_;
  if (transA == 0)
    transA_ = CblasNoTrans;
  else
    transA_=CblasTrans;
  if (transB == 0)
    transB_ = CblasNoTrans;
  else
    transB_=CblasTrans;

  gsl_blas_dgemm(transA_,transB_,alpha,A,B,beta,C);
#elif defined(LAPACK) || defined(MKL)
  int M,N,K,LDA,LDB,LDC;
  CBLAS_TRANSPOSE transA_, transB_;
  if (transA == 0){
    transA_ = CblasNoTrans;
    M = A->size1;
    K = A->size2;
    LDA = K;
  }
  else {
    transA_ = CblasTrans;
    M = A->size2;
    K = A->size1;
    LDA = M;
  }
  if (transB == 0){
    transB_ = CblasNoTrans;
    N = B->size2;
    LDB = N;
  }
  else {
    transB_ = CblasTrans;
    N = B->size1;
    LDB = K;
  }
  LDC = N;
  cblas_dgemm(CblasRowMajor, transA_, transB_,
                 M, N, K,
                 alpha, A->data, LDA, B->data, LDB,
                 beta, C->data, LDC);
#endif
  return;
}

void
linalg_dgemv  (int transA, const double alpha,
               const gsl_matrix *A, const gsl_vector *v,
               const double beta, gsl_vector *w){
#ifdef GSL
  CBLAS_TRANSPOSE_t transA_;
  if (transA == 0)
    transA_ = CblasNoTrans;
  else
    transA_=CblasTrans;

  gsl_blas_dgemv(transA_,alpha,A,v,beta,w);
#elif defined(LAPACK) || defined(MKL)
  CBLAS_TRANSPOSE transA_;
  int M,N,LDA;
  M = A->size1;
  N = A->size2;
  LDA = N;
  if (transA == 0)
    transA_ = CblasNoTrans;
  else
    transA_ = CblasTrans;

  cblas_dgemv(CblasRowMajor, transA_,
                 M, N,
                 alpha, A->data, LDA,
                 v->data, 1,
                 beta, w->data, 1);
#endif
  return;
}

void
linalg_spdgemv(const int transa, const double alpha, const  SparseMat *A, const gsl_vector *v, const double beta, gsl_vector *w){

#ifdef GSL
  CBLAS_TRANSPOSE_t transa_;
  if (transa == 0)
    transa_=CblasNoTrans;
  else
    transa_=CblasTrans;
  gsl_spblas_dgemv(transa_, alpha, A->gslmatrix, v, beta, w);
#elif defined(LAPACK)
  int N,M,K,LWORK;
  double *WORK(0);
  int *descra(0);
  M=A->nrow;
  K=A->ncol;
  LWORK=(K>M?K:M);
  WORK = new double[LWORK];
  descra = new int[9];// NIS: ftp://gams.nist.gov/pub/karin/spblas/uguide.ps.gz
  for (size_t i=0;i<9;i++) descra[i]=0;

  dcoomm(transa,M,1,K,alpha,descra,A->data,A->indx,A->indy,A->nz,v->data,v->size,beta,w->data,w->size,WORK,LWORK);

  delete[] WORK;
  delete[] descra;

#elif defined(MKL)
  char *descra(0);
  char transa_;
  int M,K;
  if (transa == 0)
    transa_ = 'N';
  else
    transa_ = 'T';
  M = A->nrow;
  K = A->ncol;
  descra = new char[6];// MKL: https://software.intel.com/en-us/node/468532#TBL2-6
  descra[0]='G';  // general matrix
  descra[1]='U';  // irrelevant (upper/lower triangular indicator)
  descra[2]='N';  // non-unitary
  descra[3]='C';  // zero-based indexing like C++

  mkl_dcoomv(&transa_,&M,&K,&alpha,descra,A->data,A->indx,A->indy,&(A->nz),v->data, &beta,w->data);
  delete[] descra;
#endif
  return;
}

void
linalg_invert(const gsl_matrix *a, gsl_matrix *ainv){

  int N;

  N = a->size1;
  if (a->size2 != N)
    throw invalid_argument("a must be a square matrix!");
#ifdef GSL
  int status;
  gsl_matrix *lu(0);
  gsl_permutation *per(0);
  status=0;

  lu = gsl_matrix_calloc(N,N);
  per = gsl_permutation_alloc(N);

  gsl_matrix_memcpy(lu,a);
  gsl_linalg_LU_decomp(lu,per,&status); // LU decomposition of input matrix
  gsl_linalg_LU_invert(lu,per,ainv);    // inverse x in xinv

  gsl_matrix_free(lu);
  gsl_permutation_free(per);

#elif defined(LAPACK) || defined(MKL)
  int INFO,LAYOUT,LDA,LWORK,MONE;
  double LWKOPT;
  int *IPIV(0);
  double *WORK(0);

  INFO = 0;
  MONE = -1;
  LAYOUT = LAPACK_ROW_MAJOR;
  LDA = N;  // row-major layout
  IPIV = new int[N];
  gsl_matrix_memcpy(ainv,a);

  INFO = LAPACKE_dgetrf(LAYOUT, N, N, ainv->data, LDA, IPIV);
  INFO = LAPACKE_dgetri_work(LAYOUT, N, ainv->data, LDA, IPIV, &LWKOPT, MONE ); // query optimal workspace size
  LWORK = int(LWKOPT);
  if (LWORK < N)  LWORK = N;
  WORK = new double[LWORK];

  INFO = LAPACKE_dgetri_work(LAYOUT, N, ainv->data, LDA, IPIV, WORK, LWORK);

  delete[] WORK;
  delete[] IPIV;
#endif
  return;
}

