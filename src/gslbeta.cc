#include <gsl/gsl_randist.h>
#include <iostream>
using namespace std;

int main()
{
  cout << gsl_ran_beta_pdf(0.5, 2, 3) << endl;
}
