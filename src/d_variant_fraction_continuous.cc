#include "setting.hh"
#include <iostream>
#include <fstream>
#include <vector>
#include <gsl/gsl_sf_gegenbauer.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_errno.h>
#include <xmmintrin.h>
using namespace std;

#define calc_gamma_i(i, n, t, beta_tilda) (i*(i+1) * beta_tilda * t / 2 / n / (log(N) + log(n)))
#define FRACTIONS 10
#define GEGEN_MAX 200
#define BETA_TILDA_MAX 10

typedef std::vector< double > V;
typedef std::vector< V > VV;
typedef double (*myfunc) (int s, int h, int q, double n_q, double t_q, double t_q_h, double beta_tilda_q, VV& gegen, V& gegen_int, double& d_denom, ofstream& f);

void write_vector(ofstream &f, V &a, int n)
{
  for (int i=0; i<n; ++i)
    f << a[i] << " ";
  f << endl;
}

void write_matrix(ofstream &f, VV &a, int m, int n)
{
  for (int i=0; i<m; ++i)
    {
      for (int j=0; j<n; ++j)
        f << a[i][j] << " ";
      f << endl;
    }
  f << endl;
}

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

double variant_fraction(int s, int h, int q, double n_q, double t_q, double t_q_h, double beta_tilda_q, VV& gegen, V& gegen_int, double& denom, ofstream& f)
{
  double x_q = ((double) s) / FRACTIONS;
  cout << "x_q = " << x_q << endl;
  double Nnq = N*n_q;
  double lNnq = log(Nnq);

  if (h == 0)
    {
      double part_acc = 0;
      for (int i=GEGEN_MAX; i>0; --i)
        {
          double g = calc_gamma_i(i, n_q, t_q, beta_tilda_q);

          part_acc += (2.0*i + 1.0) / i / (i+1.0) *
            (exp(log(g) + g + log(gsl_sf_expint_E1(g) - gsl_sf_expint_E1(g*Nnq)))  - 1.0 + exp(-g * (Nnq - 1.0))/Nnq);
        }
      denom = t_q * (1.0 - (1.0 - 1.0/Nnq)/lNnq) + 2.0 * n_q / beta_tilda_q * ((lNnq - 1.0) * (1.0 - 1.0/Nnq) + part_acc);
      f << "denom = " << denom << endl;
      
      if (s == 0)
        cout << "err: variant fraction <= 0" << endl;
      
      else if (s < FRACTIONS)
        {
          cout << "0 < x_q < 1" << endl;
          
          double acc = 0;
          for (int i=GEGEN_MAX; i>0; --i)
            {
              double g = calc_gamma_i(i, n_q, t_q, beta_tilda_q);

              double a = 
                (2.0*i + 1.0) / i / (i+1.0) * gegen[s][i] *
                (exp(log(g) + g + log(gsl_sf_expint_E1(g) - gsl_sf_expint_E1(g*Nnq)))  - 1.0 + 1.0/g + exp(-g * (Nnq - 1.0))/Nnq);
              cout << a << endl;

              acc += a;
            }
          double numerator =  2.0 * n_q / beta_tilda_q * ((1.0 - 1.0/Nnq)/x_q + 4.0*n_q/beta_tilda_q/t_q * lNnq * gegen_int[s] + 2.0*acc);

          if (numerator < 0)
            {
              cout << "h0numerator < 0!,\t 2acc = " << 2.0 * acc << endl;
              cout << "(1.0 - 1.0/Nnq)/x_q = "  << (1.0 - 1.0/Nnq)/x_q << endl;
              cout << "4.0*n_q/beta_tilda_q/t_q * log(Nnq) * gegen_int[s] = " << 4.0*n_q/beta_tilda_q/t_q * lNnq * gegen_int[s] << endl;
              cout << "4.0*n_q/beta_tilda_q/t_q * log(Nnq)  = " << 4.0*n_q/beta_tilda_q/t_q * lNnq << endl;
            }
          
          cout << "numerator = " << numerator << endl << endl;

          if (numerator < 0)
            return 0;
          
          return numerator / FRACTIONS;
        }

      else // s == FRACTIONS
        {
          double acc = 0;
          for (int i=GEGEN_MAX; i>0; --i)
            {
              double g = calc_gamma_i(i, n_q, t_q, beta_tilda_q);

              double a = (2.0*i + 1.0) / i / (i+1.0) *
                (exp(log(g) + g + log(gsl_sf_expint_E1(g) - gsl_sf_expint_E1(g*Nnq)))  - 1.0 + exp(-g * (Nnq - 1.0))/Nnq);

              if (i % 2 == 1)
                a *= -1;
              
              cout << a << endl;
              acc += a;
            }
          double numerator =  t_q * (1.0 - (1.0 - 1.0/Nnq)/lNnq) + 2.0 * n_q / beta_tilda_q * (1.0/Nnq - 1.0 + acc);
          if (numerator < 0)
            cout << "h0FRACTIONnumerator < 0!,\t acc = " << acc  << endl;
          // cout << "numerator = " << numerator << endl << endl;
          if (numerator < 0)
            return 0;
          
          return numerator;
        }
      
    }
  
  else // h > 0
    {
      denom = t_q - t_q_h;
      f << "denom = " << denom << endl;
      
      double beki = exp(t_q_h / t_q * lNnq);
      double beki_1 = exp((t_q_h/ t_q  - 1.0) * lNnq);
      
      if (0 < s && s < FRACTIONS)
        {
          double acc = 0;
          for (int i=GEGEN_MAX; i>0; --i)
            {
              double g = calc_gamma_i(i, n_q, t_q, beta_tilda_q);
              double a =
                (2.0*i + 1.0) / i / (i+1.0) * gegen[s][i] * ( exp(-g * (beki - 1.0)) - exp(-g * (Nnq - 1.0)) );
              cout << a << endl;
              acc += a;
            }
          double numerator = 4.0 * acc / N / beta_tilda_q;
          if (numerator < 0)
            cout << "numerator < 0!,\t acc = " << acc  << endl;
          // cout << "numerator = " << numerator << endl << endl;
          if (numerator < 0)
            return 0;
          
          return numerator / FRACTIONS;
        }

      else if (s == FRACTIONS)
        {
          cout << "x_q = 1" << endl;
          
          double acc = 0;
          for (int i=GEGEN_MAX; i>0; --i)
            {
              double g = calc_gamma_i(i, n_q, t_q, beta_tilda_q);
              double a =
                (2.0*i + 1.0) / i / (i+1.0) * ( exp(-g * (beki - 1.0)) - exp(-g * (Nnq - 1.0)) );
              
              if (i % 2 == 1)
                a *= -1;
              
              cout << a << endl;
              acc += a;
            }
          double numerator = t_q / lNnq * ( 1.0 - beki_1 ) + 2.0 * acc / N / beta_tilda_q;
          if (numerator < 0)
            cout << "numerator < 0!,\t acc = " << acc  << endl;
          // cout << "numerator = " << numerator << endl << endl;
          if (numerator < 0)
            return 0;
          
          return numerator;
        }

      else // s == 0
        {
          double acc = 0;
          for (int i=GEGEN_MAX; i>0; --i)
            {
              double g = calc_gamma_i(i, n_q, t_q, beta_tilda_q);
              double a =
                (2.0*i + 1.0) / i / (i+1.0) * ( exp(-g * (beki - 1.0)) - exp(-g * (Nnq - 1.0)) );
              
              cout << a << endl;
              acc += a;
            }
          double numerator = t_q - t_q_h - t_q / lNnq * ( 1.0 - beki_1 ) - 2.0 * acc / N / beta_tilda_q;
          if (numerator < 0)
            cout << "numerator < 0!,\t acc = " << acc  << endl;
          // cout << "numerator = " << numerator << endl << endl;
          if (numerator < 0)
            return 0;
          
          return numerator;
        }
    }
  return -1;
}

void variant_fraction_normalized(int h, int q, double n_q, double t_q, double t_q_h, double beta_tilda_q, VV& gegen, V& gegen_int, V& vf, double& denom, ofstream &f)
{
  double interval = 1.0 / FRACTIONS;
  if (h == 0)
    {
      for (int s=1; s<=FRACTIONS; ++s)
        {
          vf[s] = variant_fraction(s, 0, q, n_q, t_q, t_q_h, beta_tilda_q, gegen, gegen_int, denom, f);
          vf[s] /= denom;
        }
    }
  else
    {
      for (int s=0; s<=FRACTIONS; ++s)
        {
          vf[s] = variant_fraction(s, 1, q, n_q, t_q, t_q_h, beta_tilda_q, gegen, gegen_int, denom, f);
          vf[s] /= denom;
        }
    }
}

double d_t_variant_fraction(int s, int h, int q, double n_q, double t_q, double t_q_h, double beta_tilda_q, VV& gegen, V& gegen_int, double& d_denom, ofstream& f)
{
  double x_q = ((double) s) / FRACTIONS;
  cout << "x_q = " << x_q << endl;
  f << "x_q = " << x_q << endl;
  double Nnq = N*n_q;
  double lNnq = log(Nnq);

  if (h == 0)
    {
      double part_acc = 0;
      for (int i=GEGEN_MAX; i>0; --i)
        {
          double g = calc_gamma_i(i, n_q, t_q, beta_tilda_q);

          part_acc += (2.0*i + 1.0) *
            (exp(log(1.0 + g) + g + log(gsl_sf_expint_E1(g) - gsl_sf_expint_E1(g*Nnq))) - 1.0 + exp(-g * (Nnq - 1.0))/Nnq);
        }
      d_denom = 1.0 + ( 1.0/Nnq - 1.0 + part_acc ) / lNnq;
      f << "d_denom = " << d_denom << endl;
      
      if (s == 0)
        cout << "err: variant fraction <= 0" << endl;
      
      else if (s < FRACTIONS)
        {
          double acc = 0;
          for (int i=GEGEN_MAX; i>0; --i)
            {
              double g = calc_gamma_i(i, n_q, t_q, beta_tilda_q);

              double a = 
                (2.0*i + 1.0) * gegen[s][i] *
                (exp(log(1.0 + g) + g + log(gsl_sf_expint_E1(g) - gsl_sf_expint_E1(g*Nnq))) - 1.0 + 1.0/g/g + exp(-g * (Nnq - 1.0))/Nnq);
              cout << a << endl;

              acc += a;
            }
          return (-8.0 * pow(n_q/beta_tilda_q/t_q, 2.0) * lNnq * gegen_int[s] + 2.0 / lNnq * acc) / FRACTIONS;
        }

      else // s == FRACTIONS
        {
          double acc = 0;
          for (int i=GEGEN_MAX; i>0; --i)
            {
              double g = calc_gamma_i(i, n_q, t_q, beta_tilda_q);

              double a = (2.0*i + 1.0) *
                (exp(log(1.0 + g) + g + log(gsl_sf_expint_E1(g) - gsl_sf_expint_E1(g*Nnq))) - 1.0 + exp(-g * (Nnq - 1.0))/Nnq);

              if (i % 2 == 1)
                a *= -1;
              
              cout << a << endl;
              acc += a;
            }
          double numerator =  1.0 +  (1.0/Nnq - 1.0 + acc) / lNnq;
          // cout << "numerator = " << numerator << endl << endl;
          
          return numerator;
        }
      
    }
  
  else // h > 0
    {
      d_denom = 1.0;
      
      f << "d_denom = " << d_denom << endl;
      
      double beki = exp(t_q_h / t_q * lNnq);
      double beki_1 = exp((t_q_h/ t_q  - 1.0) * lNnq);
      
      if (0 < s && s < FRACTIONS)
        {
          double acc = 0;
          for (int i=GEGEN_MAX; i>0; --i)
            {
              double g = calc_gamma_i(i, n_q, t_q, beta_tilda_q);
              double a =
                (2.0*i + 1.0) * gegen[s][i] *
                ( (1.0 / lNnq / Nnq + (t_q_h / t_q - 1.0/lNnq) * beki_1) * exp(-g * (beki - 1.0)) + (1.0 - 1.0/Nnq) / lNnq * exp(-g * (Nnq - 1.0)) );
              cout << a << endl;
              acc += a;
            }
          return 2.0 * acc / FRACTIONS;
        }

      else if (s == FRACTIONS)
        {
          cout << "x_q = 1" << endl;
          
          double acc = 0;
          for (int i=GEGEN_MAX; i>0; --i)
            {
              double g = calc_gamma_i(i, n_q, t_q, beta_tilda_q);
              double a =
                (2.0*i + 1.0) *
                ( (1.0 / lNnq / Nnq + (t_q_h / t_q - 1.0/lNnq) * beki_1) * exp(-g * (beki - 1.0)) + (1.0 - 1.0/Nnq) / lNnq * exp(-g * (Nnq - 1.0)) );
              
              if (i % 2 == 1)
                a *= -1;
              
              cout << a << endl;
              acc += a;
            }
          double numerator = 1.0/lNnq + ( t_q_h / t_q - 1.0/lNnq ) * beki_1 + acc;

          // cout << "numerator = " << numerator << endl << endl;
          
          return numerator;
        }

      else // s == 0
        {
          double acc = 0;
          for (int i=GEGEN_MAX; i>0; --i)
            {
              double g = calc_gamma_i(i, n_q, t_q, beta_tilda_q);
              double a =
                (2.0*i + 1.0) *
                ( (1.0 / lNnq / Nnq + (t_q_h / t_q - 1.0/lNnq) * beki_1) * exp(-g * (beki - 1.0)) + (1.0 - 1.0/Nnq) / lNnq * exp(-g * (Nnq - 1.0)) );
              
              cout << a << endl;
              acc += a;
            }
          double numerator = 1.0 - 1.0 / lNnq  - ( t_q_h / t_q - 1.0/lNnq ) * beki_1 - acc;

          // cout << "numerator = " << numerator << endl << endl;
          
          return numerator;
        }
    }
  return -1;
}

void d_variant_fraction_normalized(myfunc d_x_variant_fraction, int h, int q, double n_q, double t_q, double t_q_h, double beta_tilda_q, VV& gegen, V& gegen_int, V& vf, V& dvf, V& ndvf, double& denom, double& d_denom, ofstream &f) // needs 
{
  double interval = 1.0 / FRACTIONS;
  double dvf_partition = 0;

  if (h == 0)
    {
      for (int s=1; s<=FRACTIONS; ++s)
        {
          dvf[s] = d_x_variant_fraction(s, 0, q, n_q, t_q, t_q_h, beta_tilda_q, gegen, gegen_int, d_denom, f);
        }

      f << "not_partitioned_d_x_variant_fraction h0:" << endl;
      write_vector(f, dvf, FRACTIONS+1);
      cout << endl;      

      for (int s=1; s<=FRACTIONS; ++s)
        {
          ndvf[s] = (dvf[s] * denom - vf[s] * d_denom) / d_denom / d_denom;
        }
    }
  else
    {
      for (int s=0; s<=FRACTIONS; ++s)
        {
          dvf[s] = d_x_variant_fraction(s, 1, q, n_q, t_q, t_q_h, beta_tilda_q, gegen, gegen_int, d_denom, f);
        }
      
      f << "not_partitioned_d_x_variant_fraction h0:" << endl;
      write_vector(f, dvf, FRACTIONS+1);
      cout << endl;      

      for (int s=0; s<=FRACTIONS; ++s)
        {
          ndvf[s] = (dvf[s] * denom - vf[s] * d_denom) / d_denom / d_denom;
        }
    }
}

int main(int argc, char **argv)
{
  _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);
  
  gsl_set_error_handler_off ();
    
  VV gegen;
  gegen.assign(FRACTIONS+1, V(GEGEN_MAX+1, 0));

  set_gegen(gegen);

  V gegen_int (FRACTIONS+1, 0);
  V gegen_int_err (FRACTIONS+1, 0);

  set_gegen_integral(gegen_int, gegen_int_err);

  V vf_h0 (FRACTIONS+1, 0);
  V vf_h1 (FRACTIONS+1, 0);

  V dvf_h0 (FRACTIONS+1, 0);
  V dvf_h1 (FRACTIONS+1, 0);

  V ndvf_h0 (FRACTIONS+1, 0);
  V ndvf_h1 (FRACTIONS+1, 0);

  ofstream f(argv[1]);

  for (int beta_tilda_disc=1; beta_tilda_disc<=BETA_TILDA_MAX; ++beta_tilda_disc)
    {
      double beta_tilda_q = ((double) beta_tilda_disc) / BETA_TILDA_MAX;

      f << "beta_tilda_q: " << beta_tilda_q << endl;
      cout << "beta_tilda_q: " << beta_tilda_q << endl;

      double denom_h0 = 0;
      double denom_h1 = 0;

      double d_denom_h0 = 0;
      double d_denom_h1 = 0;

      cout << "h = 0" << endl;
      f << "h = 0" << endl;
      variant_fraction_normalized(0, 4, 0.1, 0.2, 0.05, beta_tilda_q, gegen, gegen_int, vf_h0, denom_h0, f);
      // write_vector(f, vf_h0, FRACTIONS+1);
      d_variant_fraction_normalized(d_t_variant_fraction, 0, 4, 0.1, 0.2, 0.05, beta_tilda_q, gegen, gegen_int, vf_h0, dvf_h0, ndvf_h0, denom_h0, d_denom_h0, f);
      f << "d_t_result: ";
      write_vector(f, ndvf_h0, FRACTIONS+1);
      f << endl;

      cout << "h > 0" << endl;
      f << "h > 0" << endl;
      variant_fraction_normalized(1, 4, 0.1, 0.2, 0.05, beta_tilda_q, gegen, gegen_int, vf_h1, denom_h1, f);
      // write_vector(f, vf_h1, FRACTIONS+1);
      d_variant_fraction_normalized(d_t_variant_fraction, 0, 4, 0.1, 0.2, 0.05, beta_tilda_q, gegen, gegen_int, vf_h1, dvf_h1, ndvf_h1, denom_h1, d_denom_h1, f);
      f << "d_t_result: ";
      write_vector(f, ndvf_h1, FRACTIONS+1);

      f << endl << endl;
    }
  
  f.close();
  
  return 0;
}
