#include <iostream>
#include <fstream>
#include <vector>
#include <gsl/gsl_sf_gegenbauer.h>
using namespace std;

#define FRACTIONS 1000
#define GEGEN_MAX 1000

typedef std::vector< double > V;
typedef std::vector< V > VV;

void set_gegen(VV &gegen)
{
  for (int s=0; s<=FRACTIONS; ++s)
    {
      double frac = ((double)s) / FRACTIONS;

      for (int i=1; i<=GEGEN_MAX; ++i)
        gegen[s][i] = gsl_sf_gegenpoly_n(i-1, 1.5, 1-2*frac);
    }
}

void write_matrix(ofstream &f, VV &a, int m, int n)
{
  f << scientific;
  
  for (int i=0; i<m; ++i)
    {
      for (int j=0; j<n; ++j)
        f << a[i][j] << " ";
      f << endl;
    }
  f << endl;
}

int main()
{
  VV gegen;
  gegen.assign(FRACTIONS+1, V(GEGEN_MAX+1, 0));

  set_gegen(gegen);
  
  ofstream f;
  f.open("gegen.txt");

  write_matrix(f, gegen, FRACTIONS+1, GEGEN_MAX+1);
  
  return 0;
}
