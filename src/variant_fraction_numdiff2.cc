#include "setting.hh"
#include <iostream>
#include <fstream>
#include <vector>
#include <gsl/gsl_sf_gegenbauer.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_deriv.h>
#include <xmmintrin.h>
using namespace std;

#define calc_gamma_i(i, n, t, beta_tilda) (i*(i+1) * beta_tilda * t / 2 / n / (log(N) + log(n)))
#define FRACTIONS 10
#define GEGEN_MAX 200
#define BETA_TILDA_MAX 10
#define T_MAX 10
#define T_H_MAX 10
#define NQ_MAX 10
#define N 1000000

typedef std::vector< double > V;
typedef std::vector< V > VV;

void write_vector_vert(ofstream &f, V &a, int n)
{
  for (int i=0; i<n; ++i)
    f << ((double) i) / T_MAX << "\t" << a[i] << endl;
}

void write_vector(ofstream &f, V &a, int n)
{
  f << scientific;
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

double variant_fraction(int s, int h, int q, double n_q, double t_q, double t_q_h, double beta_tilda_q, VV& gegen, V& gegen_int)
{
  double x_q = ((double) s) / FRACTIONS;
  // cout << "x_q = " << x_q << endl;
  double Nnq = N*n_q;
  double lNnq = log(Nnq);

  if (h == 0)
    {
      double part_acc = 0;
      for (int i=GEGEN_MAX; i>0; --i)
        {
          double g = calc_gamma_i(i, n_q, t_q, beta_tilda_q);
          // cout << "i=" << i << ", g=" << g << endl;

          part_acc += (2.0*i + 1.0) / i / (i+1.0) *
            (exp(log(g) + g + log(gsl_sf_expint_E1(g) - gsl_sf_expint_E1(g*Nnq)))  - 1.0 + exp(-g * (Nnq - 1.0))/Nnq);
        }
      double partition = t_q * (1.0 - (1.0 - 1.0/Nnq)/lNnq) + 2.0 * n_q / beta_tilda_q * ((lNnq - 1.0) * (1.0 - 1.0/Nnq) + part_acc);
      // cout << "partition = " << partition << endl;
      
      if (s == 0)
        cout << "err: variant fraction <= 0" << endl;
      
      else if (s < FRACTIONS)
        {
          // cout << "0 < x_q < 1" << endl;
          
          double acc = 0;
          for (int i=GEGEN_MAX; i>0; --i)
            {
              double g = calc_gamma_i(i, n_q, t_q, beta_tilda_q);
              // cout << "i=" << i << ", g=" << g << endl;

              // cout << "(2.0*i + 1.0) / i / (i+1.0) = " << (2.0*i + 1.0) / i / (i+1.0) << endl;

              double a = 
                (2.0*i + 1.0) / i / (i+1.0) * gegen[s][i] *
                (exp(log(g) + g + log(gsl_sf_expint_E1(g) - gsl_sf_expint_E1(g*Nnq)))  - 1.0 + 1.0/g + exp(-g * (Nnq - 1.0))/Nnq);
              // cout << a << endl;

              acc += a;
            }
          double numerator =  2.0 * n_q / beta_tilda_q * ((1.0 - 1.0/Nnq)/x_q + 4.0*n_q/beta_tilda_q/t_q * lNnq * gegen_int[s] + 2.0*acc);

        //   cout << "acc1" << endl;
        //   double acc1 = 0;
        //   for (int i=GEGEN_MAX; i>0; --i)
        //     {
        //       double g = calc_gamma_i(i, n_q, t_q, beta_tilda_q);
        //       // cout << "i=" << i << ", g=" << g << endl;

        //       // cout << "(2.0*i + 1.0) / i / (i+1.0) = " << (2.0*i + 1.0) / i / (i+1.0) << endl;
              
        //       double a = 
        //         (2.0*i + 1.0) / i / (i+1.0) * gegen[s][i] *
        //         exp(log(g) + g + log(gsl_sf_expint_E1(g) - gsl_sf_expint_E1(g*Nnq)));
        //       cout << a << endl;

        //       acc1 += a;
        //     }
        //   cout << "acc1 = " << acc1 << endl;
          
        //   cout << "acc2" << endl;
        //   double acc2 = 0;
        //   for (int i=GEGEN_MAX; i>0; --i)
        //     {
        //       double g = calc_gamma_i(i, n_q, t_q, beta_tilda_q);
        //       // cout << "i=" << i << ", g=" << g << endl;

        //       // cout << "(2.0*i + 1.0) / i / (i+1.0) = " << (2.0*i + 1.0) / i / (i+1.0) << endl;
              
        //       double a = 
        //         - (2.0*i + 1.0) / i / (i+1.0) * gegen[s][i];
        //       cout << a << endl;

        //       acc2 += a;
        //     }
        //   cout << "acc2 = " << acc2 << endl;
          
        //   cout << "acc3" << endl;
        //   double acc3 = 0;
        //   for (int i=GEGEN_MAX; i>0; --i)
        //     {
        //       double g = calc_gamma_i(i, n_q, t_q, beta_tilda_q);
        //       // cout << "i=" << i << ", g=" << g << endl;

        //       // cout << "(2.0*i + 1.0) / i / (i+1.0) = " << (2.0*i + 1.0) / i / (i+1.0) << endl;
              
        //       double a = 
        //         (2.0*i + 1.0) / i / (i+1.0) * gegen[s][i] / g;
        //       cout << a << endl;

        //       acc3 += a;
        //     }
        //   cout << "acc3 = " << acc3 << endl;
          
        //   cout << "acc4" << endl;
        //   double acc4 = 0;
        //   for (int i=GEGEN_MAX; i>0; --i)
        //     {
        //       double g = calc_gamma_i(i, n_q, t_q, beta_tilda_q);
        //       // cout << "i=" << i << ", g=" << g << endl;

        //       // cout << "(2.0*i + 1.0) / i / (i+1.0) = " << (2.0*i + 1.0) / i / (i+1.0) << endl;
              
        //       double a = 
        //         (2.0*i + 1.0) / i / (i+1.0) * gegen[s][i] *
        //         exp( -g * (Nnq - 1.0) ) / N / n_q;
        //       cout << a << endl;

        //       acc4 += a;
        //     }
        //   cout << "acc4 = " << acc4 << endl;

        //   double numerator =  2.0 * n_q / beta_tilda_q * (2*acc4 + (1.0 - 1.0/Nnq)/x_q + 2*acc2 + 2*acc3 + 4.0*n_q/beta_tilda_q/t_q * log(Nnq) * gegen_int[s] + 2*acc1);

          if (numerator < 0)
            {
              // cout << "h0numerator < 0!,\t 2acc = " << 2.0 * acc << endl;
              // cout << "(1.0 - 1.0/Nnq)/x_q = "  << (1.0 - 1.0/Nnq)/x_q << endl;
              // cout << "4.0*n_q/beta_tilda_q/t_q * log(Nnq) * gegen_int[s] = " << 4.0*n_q/beta_tilda_q/t_q * lNnq * gegen_int[s] << endl;
              // cout << "4.0*n_q/beta_tilda_q/t_q * log(Nnq)  = " << 4.0*n_q/beta_tilda_q/t_q * lNnq << endl;
            }
          
          // cout << "numerator = " << numerator << endl << endl;

          if (numerator < 0)
            return 0;
          
          return numerator / FRACTIONS;

        }

      else // s == FRACTIONS
        {
          // cout << "x_q = 1" << endl;
          
          double acc = 0;
          for (int i=GEGEN_MAX; i>0; --i)
            {
              double g = calc_gamma_i(i, n_q, t_q, beta_tilda_q);
              // cout << "i=" << i << ", g=" << g << endl;

              double a = (2.0*i + 1.0) / i / (i+1.0) *
                (exp(log(g) + g + log(gsl_sf_expint_E1(g) - gsl_sf_expint_E1(g*Nnq)))  - 1.0 + exp(-g * (Nnq - 1.0))/Nnq);

              if (i % 2 == 1)
                a *= -1;
              
              // cout << a << endl;
              acc += a;
            }
          double numerator =  t_q * (1.0 - (1.0 - 1.0/Nnq)/lNnq) + 2.0 * n_q / beta_tilda_q * (1.0/Nnq - 1.0 + acc);
          if (numerator < 0)
            // cout << "h0FRACTIONnumerator < 0!,\t acc = " << acc  << endl;
          // cout << "numerator = " << numerator << endl << endl;
          if (numerator < 0)
            return 0;
          
          return numerator;

        }
      
    }
  
  else // h > 0
    {
      // double partition = t_q - t_q_h;
      // cout << "partition = " << partition << endl;

      double beki = exp(t_q_h / t_q * lNnq);
      double beki_1 = exp((t_q_h/ t_q  - 1.0) * lNnq);
      
      if (0 < s && s < FRACTIONS)
        {
          // cout << "0 < x_q < 1" << endl;
          
          double acc = 0;
          for (int i=GEGEN_MAX; i>0; --i)
            {
              double g = calc_gamma_i(i, n_q, t_q, beta_tilda_q);
              double a =
                // (2.0*i + 1.0) / i / (i+1.0) * gegen[s][i] * exp(-g * (beki - 1.0)) * ( 1.0 - exp(-g * Nnq * (1.0 - beki_1)) );
                // (2.0*i + 1.0) / i / (i+1.0) * gegen[s][i] * exp(-g * (pow(Nnq, t_q_h/t_q) - 1.0));
                // (2.0*i + 1.0) / i / (i+1.0) * gegen[s][i] * exp(-g * (pow(Nnq, t_q_h/t_q) - 1.0)) * ( 1.0 - exp(-g * Nnq * (1.0 - pow(Nnq, t_q_h/t_q - 1.0))) );
                (2.0*i + 1.0) / i / (i+1.0) * gegen[s][i] * ( exp(-g * (beki - 1.0)) - exp(-g * (Nnq - 1.0)) );
              // cout << a << endl;
              acc += a;
            }
          double numerator = 4.0 * acc / N / beta_tilda_q;
          if (numerator < 0)
            // cout << "numerator < 0!,\t acc = " << acc  << endl;
          // cout << "numerator = " << numerator << endl << endl;
          if (numerator < 0)
            return 0;
          
          return numerator / FRACTIONS;


        }

      else if (s == FRACTIONS)
        {
          // cout << "x_q = 1" << endl;
          
          double acc = 0;
          for (int i=GEGEN_MAX; i>0; --i)
          // for (int i=1; i<=GEGEN_MAX; ++i)
            {
              double g = calc_gamma_i(i, n_q, t_q, beta_tilda_q);
              double a =
                // (2.0*i + 1.0) / i / (i+1.0) * exp(-g * (beki - 1.0)) * ( 1.0 - exp(-g * Nnq * (1.0 - beki_1)) );;
                // (2.0*i + 1.0) / i / (i+1.0) * exp(-g * (pow(Nnq, t_q_h/t_q) - 1.0)) * ( 1.0 - exp(-g * Nnq * (1.0 - pow(Nnq, t_q_h/t_q - 1.0))) );
                (2.0*i + 1.0) / i / (i+1.0) * ( exp(-g * (beki - 1.0)) - exp(-g * (Nnq - 1.0)) );
              
              if (i % 2 == 1)
                a *= -1;
              
              // cout << a << endl;
              acc += a;
            }
          double numerator = t_q / lNnq * ( 1.0 - beki_1 ) + 2.0 * acc / N / beta_tilda_q;
          if (numerator < 0)
            // cout << "numerator < 0!,\t acc = " << acc  << endl;
          // cout << "numerator = " << numerator << endl << endl;
          if (numerator < 0)
            return 0;
          
          return numerator;

        }

      else // s == 0
        {
          // cout << "x_q = 0" << endl;
          
          double acc = 0;
          for (int i=GEGEN_MAX; i>0; --i)
          // for (int i=1; i<=GEGEN_MAX; ++i)
            {
              double g = calc_gamma_i(i, n_q, t_q, beta_tilda_q);
              double a =
                // (2.0*i + 1.0) / i / (i+1.0) * exp(-g * (beki - 1.0)) * ( 1.0 - exp(-g * Nnq * (1.0 - beki_1)) );
                // (2.0*i + 1.0) / i / (i+1.0) * exp(-g * (pow(Nnq, t_q_h/t_q) - 1.0)) * ( 1.0 - exp(-g * Nnq * (1.0 - pow(Nnq, t_q_h/t_q - 1.0))) );
                (2.0*i + 1.0) / i / (i+1.0) * ( exp(-g * (beki - 1.0)) - exp(-g * (Nnq - 1.0)) );
              
              // cout << a << endl;
              acc += a;
            }
          double numerator = t_q - t_q_h - t_q / lNnq * ( 1.0 - beki_1 ) - 2.0 * acc / N / beta_tilda_q;
          if (numerator < 0)
            // cout << "numerator < 0!,\t acc = " << acc  << endl;
          // cout << "numerator = " << numerator << endl << endl;
          if (numerator < 0)
            return 0;
          
          return numerator;

        }
    }
  
  return -1;
}

double variant_fraction_normalized(int s, int h, int q, double n_q, double t_q, double t_q_h, double beta_tilda_q, VV& gegen, V& gegen_int)
{
  V vf(FRACTIONS+1, 0);
  double partition = 0;

  if (h == 0)
    {
      for (int x=1; x<=FRACTIONS; ++x)
        {
          vf[x] = variant_fraction(x, 0, q, n_q, t_q, t_q_h, beta_tilda_q, gegen, gegen_int);
        }

      for (int x=1; x<=FRACTIONS; ++x)
        {
          partition += vf[x];
        }
  
      return vf[s] /= partition;
    }
  else
    {
      for (int x=0; x<=FRACTIONS; ++x)
        {
          vf[x] = variant_fraction(x, 1, q, n_q, t_q, t_q_h, beta_tilda_q, gegen, gegen_int);
        }

      for (int x=0; x<=FRACTIONS; ++x)
        {
          partition += vf[x];
        }
  
      return vf[s] /= partition;
    }
}

typedef struct _vfnumdiff_t
{
  int s;
  int h;
  int q;
  double n_q;
  double t_q_h;
  double beta_tilda_q;
  VV gegen;
  V gegen_int;
}
  vfnumdiff_t;

double func_t(double t_q, void *params)
{
  vfnumdiff_t *p = (vfnumdiff_t*) params;
  
  return variant_fraction_normalized(p->s, p->h, p->q, p->n_q, t_q, p->t_q_h, p->beta_tilda_q, p->gegen, p->gegen_int);
}

void d_t_variant_fraction_normalized(int s, int h, int q, double n_q, double t_q, double t_q_h, double beta_tilda_q, VV& gegen, V& gegen_int, double* result, double* abserr)
{
  gsl_function F;

  F.function = &func_t;

  vfnumdiff_t dt_vf;
  dt_vf.s = s;
  dt_vf.h = h;
  dt_vf.q = q;
  dt_vf.n_q = n_q;
  dt_vf.t_q_h = t_q_h;
  dt_vf.beta_tilda_q = beta_tilda_q;
  dt_vf.gegen = gegen;
  dt_vf.gegen_int = gegen_int;

  F.params = &dt_vf;
              
  gsl_deriv_backward(&F, t_q, 1e-3, result, abserr);
}

typedef struct _vfnumdiff_t_h
{
  int s;
  int h;
  int q;
  double n_q;
  double t_q;
  double beta_tilda_q;
  VV gegen;
  V gegen_int;
}
  vfnumdiff_t_h;

double func_t_h(double t_q_h, void *params)
{
  vfnumdiff_t_h *p = (vfnumdiff_t_h*) params;
  
  return variant_fraction_normalized(p->s, p->h, p->q, p->n_q, p->t_q, t_q_h, p->beta_tilda_q, p->gegen, p->gegen_int);
}

void d_t_h_variant_fraction_normalized(int s, int h, int q, double n_q, double t_q, double t_q_h, double beta_tilda_q, VV& gegen, V& gegen_int, double* result, double* abserr)
{
  gsl_function F;

  F.function = &func_t_h;

  vfnumdiff_t_h dt_h_vf;
  dt_h_vf.s = s;
  dt_h_vf.h = h;
  dt_h_vf.q = q;
  dt_h_vf.n_q = n_q;
  dt_h_vf.t_q = t_q;
  dt_h_vf.beta_tilda_q = beta_tilda_q;
  dt_h_vf.gegen = gegen;
  dt_h_vf.gegen_int = gegen_int;

  F.params = &dt_h_vf;
              
  gsl_deriv_backward(&F, t_q_h, 1e-3, result, abserr);
}

typedef struct _vfnumdiff_n
{
  int s;
  int h;
  int q;
  double t_q;
  double t_q_h;
  double beta_tilda_q;
  VV gegen;
  V gegen_int;
}
  vfnumdiff_n;

double func_n(double n_q, void *params)
{
  vfnumdiff_n *p = (vfnumdiff_n*) params;

  return variant_fraction_normalized(p->s, p->h, p->q, n_q, p->t_q, p->t_q_h, p->beta_tilda_q, p->gegen, p->gegen_int);
}

void d_n_variant_fraction_normalized(int s, int h, int q, double n_q, double t_q, double t_q_h, double beta_tilda_q, VV& gegen, V& gegen_int, double* result, double* abserr)
{
  gsl_function F;

  F.function = &func_n;

  vfnumdiff_n dn_vf;
  dn_vf.s = s;
  dn_vf.h = h;
  dn_vf.q = q;
  dn_vf.t_q = t_q;
  dn_vf.t_q_h = t_q_h;
  dn_vf.beta_tilda_q = beta_tilda_q;
  dn_vf.gegen = gegen;
  dn_vf.gegen_int = gegen_int;

  F.params = &dn_vf;

  gsl_deriv_backward(&F, n_q, 1e-5, result, abserr);
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

  ofstream f(argv[1]);

  for (int _t=1; _t <= T_MAX; ++_t)
    {
      double t_q = ((double) _t) / T_MAX;
      for (int _t_h = 0; _t_h < T_H_MAX; ++_t_h)
        {
          double t_q_h = t_q * ((double) _t_h) / T_H_MAX;
          
          string outf = "./variant_fraction_numdiff2_data/variant_fraction_numdiff2_t_h1_t" + to_string(_t) + "_t_h" + to_string(_t_h);
          ofstream ff(outf);

          for (int s = 0; s<=FRACTIONS; ++s)
            {
              double x_q = ((double) s) / FRACTIONS;
              
              double result, abserr;
              d_t_variant_fraction_normalized(s, 1, 4, 0.1, t_q, t_q_h, 0.5, gegen, gegen_int, &result, &abserr);
              ff << scientific;
              ff << x_q << "\t" << result << "\t" << abserr << endl;
            }
          ff.close();
        }
    }

  for (int _t=1; _t <= T_MAX; ++_t)
    {
      double t_q = ((double) _t) / T_MAX;
      for (int _t_h = 0; _t_h < T_H_MAX; ++_t_h)
        {
          double t_q_h = t_q * ((double) _t_h) / T_H_MAX;
          
          string outf = "./variant_fraction_numdiff2_data/variant_fraction_numdiff2_t_h_h1_t" + to_string(_t) + "_t_h" + to_string(_t_h);
          ofstream ff(outf);

          for (int s = 0; s<=FRACTIONS; ++s)
            {
              double x_q = ((double) s) / FRACTIONS;
              
              double result, abserr;
              d_t_h_variant_fraction_normalized(s, 1, 4, 0.1, t_q, t_q_h, 0.5, gegen, gegen_int, &result, &abserr);
              ff << scientific;
              ff << x_q << "\t" << result << "\t" << abserr << endl;
            }
          ff.close();
        }
    }

  for (int _t=1; _t <= T_MAX; ++_t)
    {
      double t_q = ((double) _t) / T_MAX;
      for (int _t_h = 0; _t_h < T_H_MAX; ++_t_h)
        {
          double t_q_h = t_q * ((double) _t_h) / T_H_MAX;
          
          string outf = "./variant_fraction_numdiff2_data/variant_fraction_numdiff2_n_h1_t" + to_string(_t) + "_t_h" + to_string(_t_h);
          ofstream ff(outf);

          for (int s = 0; s<=FRACTIONS; ++s)
            {
              double x_q = ((double) s) / FRACTIONS;
              
              double result, abserr;
              d_n_variant_fraction_normalized(s, 1, 4, 0.1, t_q, t_q_h, 0.5, gegen, gegen_int, &result, &abserr);
              ff << scientific;
              ff << x_q << "\t" << result << "\t" << abserr << endl;
            }
          ff.close();
        }
    }

  for (int beta_tilda_disc=1; beta_tilda_disc<=BETA_TILDA_MAX; ++beta_tilda_disc)
    {
      double beta_tilda_q = ((double) beta_tilda_disc) / BETA_TILDA_MAX;

      f << "beta_tilda_q: " << beta_tilda_q << endl;

      for (int _t = 1; _t <=T_MAX; ++_t)
        {
          double t_q = ((double) _t) / T_MAX;
          
          string outf = "./variant_fraction_numdiff2_data/variant_fraction_numdiff2_t_h0_betatilda" + to_string(beta_tilda_disc) + "_t" + to_string(_t);
          ofstream ff(outf);

          for (int s = 1; s<=FRACTIONS; ++s)
            {
              double x_q = ((double) s) / FRACTIONS;
              
              double result, abserr;
              d_t_variant_fraction_normalized(s, 0, 4, 0.1, t_q, 0.2*t_q, beta_tilda_q, gegen, gegen_int, &result, &abserr);
              ff << scientific;
              ff << x_q << "\t" << result << "\t" << abserr << endl;
            }
          ff.close();
        }

      for (int _t_h = 0; _t_h<T_MAX; ++_t_h)
        {
          double t_q_h = 0.5 * ((double) _t_h) / T_MAX;
          string outf2 = "./variant_fraction_numdiff2_data/variant_fraction_numdiff2_t_h1_betatilda" + to_string(beta_tilda_disc) + "_t" + to_string(_t_h);
          ofstream ff2(outf2);

          for (int s = 0; s<=FRACTIONS; ++s)
            {
              double x_q = ((double) s) / FRACTIONS;
              
              double result, abserr;
              d_t_variant_fraction_normalized(s, 1, 4, 0.1, 0.5, t_q_h, beta_tilda_q, gegen, gegen_int, &result, &abserr);
              ff2 << scientific;
              ff2 << x_q << "\t" << result << "\t" << abserr << endl;
            }
          ff2.close();
        }
      

      for (int _t_h = 0; _t_h < T_H_MAX; ++_t_h)
        {
          double t_q = 0.5;
          
          double t_q_h = t_q * ((double) _t_h) / T_H_MAX;
          
          string outf = "./variant_fraction_numdiff2_data/variant_fraction_numdiff2_t_h_h1_betatilda" + to_string(beta_tilda_disc) + "_t_h" + to_string(_t_h);
          ofstream ff(outf);

          for (int s = 0; s<=FRACTIONS; ++s)
            {
              double x_q = ((double) s) / FRACTIONS;
              
              double result, abserr;
              d_t_h_variant_fraction_normalized(s, 1, 4, 0.1, t_q, t_q_h, beta_tilda_q, gegen, gegen_int, &result, &abserr);
              ff << scientific;
              ff << x_q << "\t" << result << "\t" << abserr << endl;
            }
          ff.close();
        }
      
      for (int _n = 1; _n <=NQ_MAX; ++_n)
        {
          double n_q = ((double) _n) / NQ_MAX;
          
          string outf = "./variant_fraction_numdiff2_data/variant_fraction_numdiff2_n_h0_betatilda" + to_string(beta_tilda_disc) + "_n" + to_string(_n);
          ofstream ff(outf);

          for (int s = 1; s<=FRACTIONS; ++s)
            {
              double x_q = ((double) s) / FRACTIONS;
              
              double result, abserr;
              d_n_variant_fraction_normalized(s, 0, 4, n_q, 0.5, 0.1, beta_tilda_q, gegen, gegen_int, &result, &abserr);
              ff << scientific;
              ff << x_q << "\t" << result << "\t" << abserr << endl;
            }
          ff.close();

          string outf2 = "./variant_fraction_numdiff2_data/variant_fraction_numdiff2_n_h1_betatilda" + to_string(beta_tilda_disc) + "_n" + to_string(_n);
          ofstream ff2(outf2);

          for (int s = 0; s<=FRACTIONS; ++s)
            {
              double x_q = ((double) s) / FRACTIONS;
              
              double result, abserr;
              d_n_variant_fraction_normalized(s, 1, 4, n_q, 0.5, 0.1, beta_tilda_q, gegen, gegen_int, &result, &abserr);
              ff2 << scientific;
              ff2 << x_q << "\t" << result << "\t" << abserr << endl;
            }
          ff2.close();

        }
    }
  f.close();
  
  return 0;
}
