#include <gsl/gsl_sf_expint.h>
#include <iostream>
using namespace std;

#define MAX_INTER 10
#define FRACTIONS 100

int main()
{
  for (int i=1; i<=FRACTIONS; ++i)
    {
      double x =  ((double) i) * MAX_INTER / FRACTIONS;
      cout << gsl_sf_expint_E1(x) << endl;
    }
  return 0;
}
