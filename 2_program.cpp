#include <armadillo>
#include <cmath>
#include <cstdlib>


using namespace arma;
using namespace std;

// Functions for Jacobi Algorithm
////////////////////////////////////////////////////
int findmax_offdiag(mat &A, int &k, int &l, double &epsilon, int n)
{
  k = 1;
  l = 0;
  double max_val = abs(A(k,l));
  
  for(int i = 1; i<n;i++)
    {
      for(int j = 0; j<i ;j++)
        {
          if (abs(A(j,i))>max_val)
             {
               max_val = abs(A(j,i));
               k = j;
               l = i;
             }
        }
    }
  if (max_val < epsilon)
    {
      return 1;
    }
  return 0;
}


void rotate(mat &A,int k, int l, int N)
{
  double cos, sin;
  if (A(k,l) != 0.0)
    {
    double tau = (A(l,l) - A(k,k))/(2*A(k,l));
    double tan;

    if (tau > 0)
      {
        tan = 1.0/(tau + sqrt(tau*tau +1));
      }
    else
      {
        tan = 1.0/(tau - sqrt(tau*tau +1));
      }
    
    cos = 1.0/(sqrt(1+tan*tan)); //Hvorfor + og ikke -?
    sin = cos*tan;
      
    }
  else
    {
      cos = 1.0;
      sin = 0.0;
    }
  
  double ik,il;

  for(int i = 0; i<N; i++)
    {
      if (i!=k && i != l)
        {
          ik = A(i,k);
          il = A(i,l);            
          A(i,k) = ik*cos - il*sin;
          A(k,i) = A(i,k);
          A(i,l) = il*cos + ik*sin;
          A(l,i) = A(i,l);
        }
    }
  
  
  double ll = A(l,l);
  double kk = A(k,k);
  double kl = A(k,l);
  
  
  A(k,k) = kk*cos*cos - 2*kl*cos*sin + ll*sin*sin;
  A(l,l) = ll*cos*cos + 2*kl*cos*sin + kk*sin*sin;
  A(k,l) = 0.0;
  A(l,k) = 0.0;  
  
}

void Jacobi_alg( mat &A, int N, double epsilon, int max_iter)
{
  int iter=0;
  int k,l,done;
  done = findmax_offdiag(A,k,l,epsilon,N); 
  while (!done && iter <= max_iter)
    {
      rotate(A,k,l,N);
      done = findmax_offdiag(A,k,l,epsilon,N);
      iter++;
    }
  if (iter == max_iter)
    {cout << "Warning: Maximum number of iterations reached"<<endl;}
}
//////////////////////////////////////////////////////////


vec HO_potential(vec rho){return rho%rho;}


void fill_matrix(mat &A, vec pot, double h)
{
  A.diag() = (2/(h*h) + pot);
  A.diag(1).fill(-1/(h*h));
  A.diag(-1).fill(-1/(h*h));
}

int main()
{
  
  int n_last_eigenvalue = 3;
  int N = 300; // N = n_{step}-1
  double rho_max = sqrt(4*n_last_eigenvalue+3)+1;
  double h = rho_max/(N+1);//skal det vaere n eller n+1 eller 2?
  vec rho(N);
  
  for(int i = 0; i<N; i++)
    {
      rho(i) = (i+1)*h;
    }
  
  mat A = zeros<mat>(N,N);
  vec pot = HO_potential(rho);
  fill_matrix(A, pot, h);
  // cout<< rho<<endl;
  // cout << pot<<endl;
  // cout << A;
  
  
  int max_iter = N*N*N;
  double epsilon = 1.0e-8;
  Jacobi_alg(A, N, epsilon, max_iter);
  
  vec eig_vals = diagvec(A);
  
  eig_vals = sort(eig_vals);
  cout<<eig_vals(span(0,n_last_eigenvalue -1)); //Noe gaar til helvete med verdiene
  
  // mat B = zeros<mat>(N,N);
  // fill_matrix(B, pot, h);
  // vec arm_eig_vals = eig_sym(A);
  // arm_eig_vals = sort(arm_eig_vals);
  // cout<<endl<<arm_eig_vals(span(0,3));
    
  return 0;
}
