#include <armadillo>


using namespace arma;


vec HO_potential(vec rho){return rho%rho;}


void fill_matrix(mat &A, vec pot, double h)
{
  A.diag() = (2/(h*h) + pot);
  A.diag(1).fill(-1/(h*h));
  A.diag(-1).fill(-1/(h*h));
}
