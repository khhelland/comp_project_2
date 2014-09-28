#pragma once
#include <armadillo>
using namespace arma;
vec HO_potential(vec rho);
void fill_matrix(mat &A, vec pot, double h);
vec two_particle_potential(vec rho, double omega);
