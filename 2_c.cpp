#include <armadillo>
#include <cmath>
#include <cstdlib>
#include "jacobi.h"
#include "matrixfunctions.h"

using namespace arma;
using namespace std;




int main()
{
  
  int N = 100; // N = n_{step}-1
  double omega = 0.25;
  double rho_max = 10;
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


  // mat B = zeros<mat>(N,N);
  // fill_matrix(B, pot, h);
  // vec arm_eig_vals = eig_sym(A);
  // arm_eig_vals = sort(arm_eig_vals);
  // cout<<endl<<arm_eig_vals(span(0,n_last_eigenvalue-1));
    
  return 0;
}
