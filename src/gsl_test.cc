#include <iostream>
#include <gsl/gsl_sf_gegenbauer.h>
#include <gsl/gsl_randist.h>
using namespace std;

int main()
{
  cout << gsl_sf_gegenpoly_n(2, 1.5, 1) << endl;
  cout << gsl_ran_binomial_pdf(1, 0.3, 2) << endl;
  cout << gsl_ran_beta_pdf(0.4, 2, 3) << endl;
  return 0;
}
