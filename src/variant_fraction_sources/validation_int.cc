#include "setting.hh"
#include <iostream>
#include <fstream>
#include <vector>
#include <gsl/gsl_sf_gegenbauer.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_errno.h>
using namespace std;

#define calc_gamma_i(i, n, t, beta_tilda) (i*(i+1) * beta_tilda * t / 2 / n / log(N * n))
#define FRACTIONS 10
#define GEGEN_MAX 200

typedef std::vector< double > V;
typedef std::vector< V > VV;

void set_gegen(VV &gegen)
{
  for (int s=1; s<=FRACTIONS; ++s)
    {
      double frac = ((double)s) / FRACTIONS;
      
      for (int i=1; i<=GEGEN_MAX; ++i)
        gegen[s][i] = gsl_sf_gegenpoly_n(i-1, 1.5, 1.0-2.0*frac);
    }
}

double f (double s, void * params) {
  double z = *(double *) params;
  double f = (1.0 - s) * log(s) / pow(1.0 - 2.0*z*s + s*s, 1.5);
  return f;
}

void set_gegen_integral(V &gegen_int, V &gegen_int_err)
{
  gsl_integration_workspace *w = gsl_integration_workspace_alloc (1000);

  double result, error;

  gsl_function F;
  F.function = &f;

  for (int s=1; s<=FRACTIONS; ++s)
    {
      double z = 1.0 - 2.0 * ((double) s) / FRACTIONS;
      
      F.params = &z;
  
      gsl_integration_qags (&F, 0, 1.0, 0, 1e-7, 1000, w, &gegen_int[s], &gegen_int_err[s]);

    }
  gsl_integration_workspace_free (w);
}



void num_err(int s, double n_q, double t_q, double beta_tilda_q, VV& gegen, V& gegen_int)
{
  double x_q = ((double) s) / FRACTIONS;
  cout << "x_q = " << x_q << endl;
  
  double acc = 0;
  // for (int i=1; i<=GEGEN_MAX; ++i)
  for (int i=GEGEN_MAX; i>0; --i)
    {
      double a = (2.0*i + 1.0) / i / (i+1.0) * gegen[s][i];
      cerr << a << endl;
      acc += a;
    }
  cerr << endl;
  cout << "-2.0*acc: " << -2.0*acc << endl;
  
  double err = 1.0/x_q - 2.0 * acc;
  double relerr = err / 2.0 / acc;
  
  cout << "err = " << err << endl;
  cout << "relerr = " << relerr << endl << endl;
}

int main()
{
  gsl_set_error_handler_off ();
    
  VV gegen;
  gegen.assign(FRACTIONS+1, V(GEGEN_MAX+1, 0));

  set_gegen(gegen);

  V gegen_int (FRACTIONS+1, 0);
  V gegen_int_err (FRACTIONS+1, 0);

  set_gegen_integral(gegen_int, gegen_int_err);

  for (int s=0; s<=FRACTIONS; ++s)
    {
      num_err(s, 0.1, 0.5, 0.05, gegen, gegen_int);
    }
  return 0;
}
