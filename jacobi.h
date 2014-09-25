#pragma once
#include <armadillo>
using namespace arma;

int Jacobi_alg( mat &A, mat &R, int N, double epsilon, int max_iter);
