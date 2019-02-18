#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <gsl/gsl_integration.h>
using namespace std;

#define FRACTIONS 10

typedef std::vector< double > V;

double f (double s, void * params) {
  double z = *(double *) params;
  double f = (1 - s) * log(s) / pow(1 - 2*z*s + s*s, 1.5);
  return f;
}

void set_gegen_integral(V &gegen_int, V &gegen_int_err)
{
  gsl_integration_workspace *w = gsl_integration_workspace_alloc (1000);

  double result, error;

  gsl_function F;
  F.function = &f;

  // double pts[2] = 
  //   {
  //     0, 1
  //   }
  // ;
  
  for (int s=1; s<=FRACTIONS; ++s)
    {
      double z = 1 - 2 * ((double) s) / FRACTIONS;
      
      F.params = &z;
  
      gsl_integration_qags (&F, 0, 1, 0, 1e-7, 1000,
                            w, &gegen_int[s], &gegen_int_err[s]);

      // printf ("s = %d\n", s);
      // printf ("result          = % .18f\n", result);
      // printf ("estimated error = % .18f\n", error);
      // printf ("intervals =  %d\n\n", (int)w->size);

    }
  gsl_integration_workspace_free (w);
}

void write_vector(ofstream &f, V &a, int n)
{
  for (int i=0; i<n; ++i)
    f << a[i] << " ";
  f << endl;
}

int main()
{
  // double* gegen_int = calloc(FRACTIONS+1, sizeof(double));
  // double* gegen_int_err = calloc(FRACTIONS+1, sizeof(double));
  
  V gegen_int (FRACTIONS+1, 0);
  V gegen_int_err (FRACTIONS+1, 0);

  set_gegen_integral(gegen_int, gegen_int_err);

  ofstream f ("numerical_integration.txt");

  write_vector(f, gegen_int, FRACTIONS+1);
  write_vector(f, gegen_int_err, FRACTIONS+1);

  return 0;
}
