//****************************************************************************
// utils_gsl.h - Defines a few useful functions.
// Date: 2017-05-03
// Created by: Guillaume Le Treut
//****************************************************************************

#ifndef UTILS_GSL_H
#define UTILS_GSL_H

#include <stdio.h>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <fstream>
#include <sstream>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_rng.h>

#include "sparse_matrix.h"

void symmatrix_vector(const gsl_matrix *m, gsl_vector *v);
void vector_symmatrix(const gsl_vector *v,gsl_matrix *m);
template <class T> void print_array(const T *a, int n);
void print_vector(const gsl_vector *v);
void print_matrix(std::ostream &mystream, const gsl_matrix *m);
void print_matrix(std::string myfile, const gsl_matrix *m);
void print_matrix(const gsl_matrix *m);
void print_matrix(const gsl_spmatrix *m);
void print_matrix(SparseMat *m);
void generate_random_matrix(gsl_matrix *m, gsl_rng *rng);
void generate_random_matrix2(gsl_matrix *m, double eps, gsl_rng *rng);
void load_matrix(std::istream &mystream, gsl_matrix *m);

#endif
