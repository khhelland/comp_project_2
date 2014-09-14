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
           if (A(j,i)>max_val)
             {
               max_val = abs(A(j,i));
               k = j;
               l = i;
             }
        }
    }
  if (max_val < epsilon)
    {
      return 0;
    }
  return 1;
}

void rotate(mat &A,int k, int l, int N)
{
  double cos, sin;
  if (A(k,l) != 0.0)
    {
    double tau = (A(l,l) - A(k,l))/(2*A(k,l));
    double tan;

    if (tau > 0)
      {
        tan = 1/(tau + sqrt(pow(tau,2) +1));
      }
    else
      {
        tan = 1/(tau - sqrt(tau*tau +1));
      }
    
    cos = 1/(sqrt(1+tan*tan)); //Hvorfor + og ikke -?
    sin = cos*tan;
      
    }
  else
    {
      cos = 1;
      sin = 0;
    }
  
  
  for(int i = 0; i<N; i++)
    {
      if (i==k && i == l)
        {
          A(i,k) = A(i,k)*cos - A(i,l)*sin;
          A(k,i) = A(i,k);
          A(i,l) = A(i,l)*cos - A(i,k)*sin;
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

  



int main()
{
  int N = 4;
  mat A(N,N);
  A.fill(2);
  int k,l, done;
  double epsilon = 1.0e-8;
  while (findmax_offdiag(A,k,l,epsilon,N))
    {rotate(A,k,l,N);}
                      
  cout<< A;
    
  return 0;
}
