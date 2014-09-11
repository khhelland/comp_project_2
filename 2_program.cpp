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
  if (max_val*max_val < epsilon)
  {
    return 1;
  }
  return 0;
}

int main()
{
  int N = 3;
  mat A = randu<mat>(N,N);
  int k,l, done;
  double epsilon = 0;
  done = findmax_offdiag(A,k,l,epsilon,N);
  cout<<A(k,l)<<endl;
  cout<<A;
  return 0;
}
