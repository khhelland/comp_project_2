#include <armadillo>
#include <cmath>
#include <cstdlib>


using namespace arma;
using namespace std;


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
    
    cos = 1.0/(sqrt(1+tan*tan)); 
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

int Jacobi_alg( mat &A, int N, double epsilon, int max_iter)
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
  return iter;
}

