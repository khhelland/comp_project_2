#include <armadillo>
#include <cmath>
#include <cstdlib>
#include "jacobi.h"
#include "matrixfunctions.h"

using namespace arma;
using namespace std;

int main()
{    
  int n_eigenvalues = 3;
  
  int N = 200; // N = n_{step}-1
  double rho_max = 4.8;
  double h = rho_max/(N+1);
  
  //create rho-vector
  vec rho(N);
  for(int i = 0; i<N; i++)
    {rho(i) = (i+1)*h;}
  
  //create matrix
  mat A = zeros<mat>(N,N);
  mat R(N,N);
  vec pot = HO_potential(rho);
  fill_matrix(A, pot, h);
  
  //maximum number of iterations
  int max_iter = N*N*N;
  //threshold for considering values zero:
  double epsilon = 1.0e-8;
  
  //Perform diagonalization
  int iter = Jacobi_alg(A, R, N, epsilon, max_iter);
  
  //print n smallest eigenvalues
  vec eig_vals = diagvec(A);
  eig_vals = sort(eig_vals);
  cout<<eig_vals(span(0,n_eigenvalues -1)); 
  
  cout<<"Number of iterations: "<<iter<<endl;


  return 0;
}
