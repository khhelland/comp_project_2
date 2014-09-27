#include <armadillo>
#include <cmath>
#include <cstdlib>
#include "jacobi.h"
#include "matrixfunctions.h"

using namespace arma;
using namespace std;




int main(int argc, char *argv[])
{
  
  int N = 300; // N = n_{step}-1
  double omega = atof(argv[1]);
  double rho_max = atof(argv[2]);
  double h = rho_max/(N+1);
  vec rho(N);
  
  for(int i = 0; i<N; i++)
    {
      rho(i) = (i+1)*h;
    }
  
  mat A = zeros<mat>(N,N);
  mat R = zeros<mat>(N,N);
  vec pot = two_independent_particle_potential(rho, omega);
  fill_matrix(A, pot, h);
  
  
  int max_iter = N*N*N;
  double epsilon = 1.0e-8;
  int iter = Jacobi_alg(A, R, N, epsilon, max_iter);
  
  vec eig_vals = diagvec(A);
  uvec indices = sort_index(eig_vals);

  
  cout<<eig_vals(indices(0))<<endl; 
  
  cout<<"Number of iterations: "<<iter<<endl;
  vec ground =R(span(),indices(0));
  ground /= sqrt(h*sum(ground%ground));
  ground.save("ground.dat",raw_ascii);
      
  return 0;
}
