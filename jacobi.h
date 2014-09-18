#pragma once
#include <armadillo>
using namespace arma;

void rotate(mat &A,int k, int l, int N);
int findmax_offdiag(mat &A, int &k, int &l, double &epsilon, int n);
int Jacobi_alg( mat &A, int N, double epsilon, int max_iter);
