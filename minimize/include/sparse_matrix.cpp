#include "sparse_matrix.h"
#include <ctime>

#ifdef GSL
SparseMat_struct::SparseMat_struct(int N, int P) :
  nrow(N),
  ncol(P)
{
  gslmatrix=gsl_spmatrix_alloc(N,P);
}

SparseMat_struct::~SparseMat_struct()
{
  gsl_spmatrix_free(gslmatrix);
}

double
SparseMat_struct::get(int i, int j){
  return gsl_spmatrix_get(gslmatrix,i,j);
}

void
SparseMat_struct::set(int i, int j, double val){
  gsl_spmatrix_set(gslmatrix,i,j,val);
  return;
}

void
SparseMat_struct::set_zero(){
  gsl_spmatrix_set_zero(gslmatrix);
  return;
}

#else

SparseMat_struct::SparseMat_struct(int N, int P) :
  nrow(N),
  ncol(P)
{
  init();
}

/* destructor */
SparseMat_struct::~SparseMat_struct()
{
  delete[] indx;
  delete[] indy;
  delete[] data;
}

/* methods */
void
SparseMat_struct::init(){

  nz=0;
  nzmax=10*std::max(nrow,ncol);    // guess for sparcity density
  if (nzmax == 0) nzmax=1;
  tol=2.0;
  indx=new int[nzmax];
  indy=new int[nzmax];
  data=new double[nzmax];
  return;
}

void
SparseMat_struct::set_zero(){
  delete[] indx;
  delete[] indy;
  delete[] data;
  init();
  return;
}

void
SparseMat_struct::extend_storage(int newnzmax){
  int *newindx(0), *newindy(0);
  double *newdata(0);

  if (newnzmax < nzmax)
    return;

  /* create new array */
  newindx = new int[newnzmax];
  newindy = new int[newnzmax];
  newdata = new double[newnzmax];

  /* copy old arrays */
  for (size_t i=0; i<nzmax; i++){
    newindx[i]=indx[i];
    newindy[i]=indy[i];
    newdata[i]=data[i];
  }

  /* delete old arrays */
  delete[] indx;
  delete[] indy;
  delete[] data;

  /* keep newly created pointer */
#ifdef DEBUG
  printf("re-allocating: nzmax=%d  nzmaxnew=%d  tol=%.2f\n",nzmax,newnzmax,tol);
#endif
  indx=newindx;
  indy=newindy;
  data=newdata;
  nzmax=newnzmax;

  return;
}

double
SparseMat_struct::get(int i, int j){
  for (size_t k=0;k<nz;k++){
    if ( (indx[k] == i) && (indy[k] == j) )
      return data[k];
  }

  return 0.0;
}

void
SparseMat_struct::set(int i, int j, double val){

  /* case where i and j are already in the matrix */
//  clock_t tinit,tter;
//  double tdiff;
//  tinit = clock();
//  for (size_t k=0;k<nz;k++){
//    if ( (indx[k] == i) && (indy[k] == j) ){
//      data[k]=val;
//      return;
//    }
//  }
//  tter = clock();
//  tdiff = (tter-tinit);
//  printf("Lookup time:  %.6e\n", tdiff/CLOCKS_PER_SEC);

  /* case where it is a new element */
  if (nz == nzmax)
    extend_storage(nzmax*tol);

  indx[nz]=i;
  indy[nz]=j;
  data[nz]=val;
  nz++;

  return;
}
#endif
