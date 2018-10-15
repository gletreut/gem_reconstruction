//****************************************************************************
// utils_gsl.cpp - Implementation of the definitions in utils_gsl.h.
// Date: 2017-05-03
// Created by: Guillaume Le Treut
//****************************************************************************
#include "utils_gsl.h"

using namespace std;

void
symmatrix_vector(const gsl_matrix *m, gsl_vector *v){
  int n1,n2,n,mu;
  n1=m->size1;
  n2=m->size2;
  n=v->size;

  if (n1 != n2)
    throw invalid_argument("m must be square");
  if ( n != (n1*(n1+1)/2) )
    throw invalid_argument("mismatch between matrix and vector");

  mu=0;
  for (int i=0;i<n1;i++)
    for (int j=i;j<n1;j++){
      gsl_vector_set(v,mu,gsl_matrix_get(m,i,j));
      mu++;
      }

  return;
}

void
vector_symmatrix(const gsl_vector *v,gsl_matrix *m){
  int n1,n2,n,mu;
  n1=m->size1;
  n2=m->size2;
  n=v->size;

  if (n1 != n2)
    throw invalid_argument("m must be square");
  if ( n != (n1*(n1+1)/2) )
    throw invalid_argument("mismatch between matrix and vector");

  mu=0;
  for (int i=0;i<n1;i++)
    for (int j=i;j<n1;j++){
      gsl_matrix_set(m,i,j,gsl_vector_get(v,mu));
      gsl_matrix_set(m,j,i,gsl_vector_get(v,mu));
      mu++;
      }

  return;
}

template <class T> void
print_array(const T *a, int n){
  cout << scientific << setprecision(4) << left;
  for (int i=0;i<n;i++){
    cout << setw(12) << a[i];
  }
  cout << endl;
  return;
}

void
print_vector(const gsl_vector *v){
  int n;
  n=v->size;

  for (int i=0;i<n;i++){
    printf("%+-12.4e",gsl_vector_get(v,i));
//    printf("%+-20.8e",gsl_vector_get(v,i));
  }
  printf("\n");
  return;
}

void
print_matrix(const gsl_matrix *m){
  int n1,n2;
  n1=m->size1;
  n2=m->size2;

  for (int i=0;i<n1;i++){
    for (int j=0;j<n2;j++){
     printf("%+-12.4e",gsl_matrix_get(m,i,j));
    }
    printf("\n");
  }
//  printf("\n");
  return;
}

void
print_matrix(ostream &mystream, const gsl_matrix *m){
  int n1,n2;
  n1=m->size1;
  n2=m->size2;

  mystream << left << showpos << scientific << setprecision(8);

  for (int i=0;i<n1;i++){
    for (int j=0;j<n2;j++){
      mystream << setw(20) << i;
      mystream << setw(20) << j;
      mystream << setw(20) << gsl_matrix_get(m,i,j);
      mystream << endl;
    }
//    mystream << endl;
  }
  return;
}

void
print_matrix(string myfileout, const gsl_matrix *m){

  ofstream fout;
  fout.open(myfileout.c_str());
  print_matrix(fout, m);
  fout.close();
  return;
}

void
print_matrix(const gsl_spmatrix *m){
  int n1,n2;
  n1=m->size1;
  n2=m->size2;

  for (int i=0;i<n1;i++){
    for (int j=0;j<n2;j++){
     printf("%+-12.4e",gsl_spmatrix_get(m,i,j));
    }
    printf("\n");
  }
  printf("\n");
  return;
}

void
print_matrix(SparseMat *m){
  int n1,n2;
  n1=m->nrow;
  n2=m->ncol;

  for (int i=0;i<n1;i++){
    for (int j=0;j<n2;j++){
     printf("%+-12.4e",m->get(i,j));
    }
    printf("\n");
  }
  printf("\n");
  return;
}

void
generate_random_matrix(gsl_matrix *m, gsl_rng *rng){
  /* generate a random matrix with positive diagonal */

  double r;
  int n1,n2;
  n1=m->size1;
  n2=m->size2;

  for (int i=0;i<n1;i++){
    r=gsl_rng_uniform(rng);
    gsl_matrix_set(m,i,i,r);
    for (int j=i+1;j<n2;j++){
      r=2*gsl_rng_uniform(rng)-1;
      gsl_matrix_set(m,i,j,r);
      gsl_matrix_set(m,j,i,r);
    }
  }
  return;
}

void
generate_random_matrix2(gsl_matrix *m, double eps, gsl_rng *rng){
  /*
   * generate a symmetric random matrix.
   * entries are uniform random variables with std eps.
   */

  double r, umax;
  int n1,n2;
  n1=m->size1;
  n2=m->size2;

  umax=sqrt(3.0)*eps;  // random variable in range -umax,+umax

  for (int i=0;i<n1;i++){
    for (int j=i;j<n2;j++){
      r=gsl_rng_uniform(rng);
      r=umax*(2*r-1.0);
      gsl_matrix_set(m,i,j,r);
      gsl_matrix_set(m,j,i,r);
    }
  }
  return;
}

void
load_matrix(std::istream &mystream, gsl_matrix *m){
	std::string line;
	std::stringstream convert;
	int i,j,n,p;
	double val;

  n=m->size1;
  p=m->size2;

	while (getline(mystream,line)){
		convert.clear();
    convert.str("");
		convert.str(line);

		if ( (convert >> i) && (convert >> j) && (convert >> val) ){
			if ( (i<n) && (j<p) && !(i<0) && !(j<0) )
				gsl_matrix_set(m,i,j,val);
		}
	}

#if defined(DEBUG_UTIL)
	std::cout << "MAT WITH SIZE " << n << " x " << p << " IMPORTED" << std::endl;
#endif
  return;
}

