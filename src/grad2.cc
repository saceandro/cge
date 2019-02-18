#include "em_io.hh"
#include <gsl/gsl_sf_gegenbauer.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_deriv.h>
#include <xmmintrin.h>
using namespace std;

#define BETA_TILDA_MAX 10
#define T_MAX 10
#define T_H_MAX 10
#define NQ_MAX 10

#define calc_gamma_i(i, n, t, beta_tilda) (i*(i+1) * beta_tilda * t / 2 / n / (log(CELL_MAX) + log(n)))

VV gegen;
V gegen_int;
V gegen_int_err;

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

double variant_fraction(int s, int h, int q, double n_q, double t_q, double t_q_h, double beta_tilda_q)
{
  double x_q = ((double) s) / FRACTIONS;
  // cout << "x_q = " << x_q << endl;
  double Nnq = CELL_MAX*n_q;
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
          double numerator = 4.0 * acc / CELL_MAX / beta_tilda_q;
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
          double numerator = t_q / lNnq * ( 1.0 - beki_1 ) + 2.0 * acc / CELL_MAX / beta_tilda_q;
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
          double numerator = t_q - t_q_h - t_q / lNnq * ( 1.0 - beki_1 ) - 2.0 * acc / CELL_MAX / beta_tilda_q;
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

double variant_fraction_normalized(int s, int h, int q, double n_q, double t_q, double t_q_h, double beta_tilda_q)
{
  V vf(FRACTIONS+1, 0);
  double partition = 0;

  if (h == 0)
    {
      for (int x=1; x<=FRACTIONS; ++x)
        {
          vf[x] = variant_fraction(x, 0, q, n_q, t_q, t_q_h, beta_tilda_q);
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
          vf[x] = variant_fraction(x, 1, q, n_q, t_q, t_q_h, beta_tilda_q);
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
}
  vfnumdiff_t;

double func_t(double t_q, void *params)
{
  vfnumdiff_t *p = (vfnumdiff_t*) params;
  
  return variant_fraction_normalized(p->s, p->h, p->q, p->n_q, t_q, p->t_q_h, p->beta_tilda_q);
}

void d_t_variant_fraction_normalized(int s, int h, int q, double n_q, double t_q, double t_q_h, double beta_tilda_q, double* result, double* abserr)
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
}
  vfnumdiff_t_h;

double func_t_h(double t_q_h, void *params)
{
  vfnumdiff_t_h *p = (vfnumdiff_t_h*) params;
  
  return variant_fraction_normalized(p->s, p->h, p->q, p->n_q, p->t_q, t_q_h, p->beta_tilda_q);
}

void d_t_h_variant_fraction_normalized(int s, int h, int q, double n_q, double t_q, double t_q_h, double beta_tilda_q, double* result, double* abserr)
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
}
  vfnumdiff_n;

double func_n(double n_q, void *params)
{
  vfnumdiff_n *p = (vfnumdiff_n*) params;

  return variant_fraction_normalized(p->s, p->h, p->q, n_q, p->t_q, p->t_q_h, p->beta_tilda_q);
}

void d_n_variant_fraction_normalized(int s, int h, int q, double n_q, double t_q, double t_q_h, double beta_tilda_q, double* result, double* abserr)
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

  F.params = &dn_vf;

  gsl_deriv_backward(&F, n_q, 1e-5, result, abserr);
}


void init_state(state* st)
{
  st->l.assign(MAX_SUBTYPE + 1, 0);
  st->r.assign(MAX_SUBTYPE + 1, 0);
}

void init_state_all(state* st)
{
  st->k = 1;
  st->l.assign(MAX_SUBTYPE + 1, 5);
  st->r.assign(MAX_SUBTYPE + 1, 2);
  st->q = 1;
  st->h = 2;
  st->s = 0.5;
  st->resp = 0.2;
}

void init_all_states(states* sts)
{
  for (states::iterator st = sts->begin(); st != sts->end(); ++st)
    {
      init_state_all(&(*st));
    }
}

void init_params(params* pa)
{
  pa->u.assign(MAX_SUBTYPE + 2, 0.5);
  pa->u[1] = 1.0;
  // for (int i=2; i<=MAX_SUBTYPE; ++i)
  //   pa->u[i] = 0.5;

  pa->beta.assign(NONLEAF + 1, std::vector<double>(MAX_CHILD, 0.5));

  pa->xi.assign(MAX_SUBTYPE, 1.0 / MAX_SUBTYPE);

  pa->omega.assign(MAX_H, 1.0 / MAX_SUBTYPE);
  
  pa->beta_tilda.assign(MAX_SUBTYPE, 0.5);

  pa->pi.assign(MAX_COPY + 1, 1.0/MAX_COPY);

  pa->kappa.assign(MAX_COPY + 1, std::vector<double>(MAX_COPY + 1, 0));
  for (int c=1; c<=MAX_COPY; ++c)
    {
      for (int d=1; d<=c; ++d)
        pa->kappa[c][d] = 1.0/c;
    }
  
  return;
}

void init_hyperparams(hyperparams* hpa)
{
  hpa->au = 2;
  hpa->bu = 2;
  hpa->abeta = 2;
  hpa->bbeta = 3;
  hpa->epsilon = 0.05;
  
  return;
}

void write_params(std::ofstream& f, params pa)
{
  f << "u" << endl;
  for (int i=2; i<=MAX_SUBTYPE; ++i)
    f << pa.u[i] << "\t";
  f << endl;
  
  f << endl << "beta" << endl;
  f << pa.beta[0][0] << endl;
  for (int i=1; i<=NONLEAF; ++i)
    {
      for (int j=0; j<MAX_CHILD; ++j)
        f << pa.beta[i][j] << "\t";
      f << endl;
    }
  
  f << endl << "pi" << endl;
  for (int c=1; c<=MAX_COPY; ++c)
    f << pa.pi[c] << "\t";
  f << endl;

  f << endl << "kappa" << endl;
  for (int c=1; c<=MAX_COPY; ++c)
    {
      for (int d=1; d<=c; ++d)
        f << pa.kappa[c][d] << "\t";
      f << endl;
    }
  
  f << endl << "xi" << endl;
  for (int i=1; i<=MAX_SUBTYPE; ++i)
    f << pa.xi[i] << "\t";
  f << endl;

  f << endl << "omega" << endl;
  for (int h=0; h<MAX_H; ++h)
    f << pa.omega[h] << "\t";
  f << endl;
  
  f << endl << "beta_tilda" << endl;
  for (int i=1; i<=MAX_SUBTYPE; ++i)
    f << pa.beta_tilda[i] << "\t";
  f << endl << endl;
  
  return;
}

void delete_reads(READS* res)
{
  for (READS::iterator it=res->begin(); it!=res->end(); ++it)
    delete *it;
  delete res;
}

bool check_converge(params pa, params prev_pa,  double eps)
{
  for (int i=2; i<=MAX_SUBTYPE; ++i)
    {
      if (fabs(pa.u[i] - prev_pa.u[i]) > eps)
        {
          cerr << "pa yet converged" << endl;
          return false;
        }
    }

  if (fabs(pa.beta[0][0] - prev_pa.beta[0][0]) > eps)
    {
      cerr << "beta yet converged" << endl;
      return false;      
    }
  
  for (int i=1; i<=NONLEAF; ++i)
    for (int j=0; j<MAX_CHILD; ++j)
      if (fabs(pa.beta[i][j] - prev_pa.beta[i][j]) > eps)
        {
          cerr << "beta yet converged" << endl;
          return false;
        }

  for (int c=1; c<=MAX_COPY; ++c)
    if (fabs(pa.pi[c] - prev_pa.pi[c]) > eps)
      {
        cerr << "pi yet converged" << endl;
        return false;
      }
  

  for (int c=1; c<=MAX_COPY; ++c)
    for (int d=1; d<=c; ++d)
      if (fabs(pa.kappa[c][d] - prev_pa.kappa[c][d]) > eps)
        {
          cerr << "kappa yet converged" << endl;        
          return false;
        }
  
  for(int i=1; i<=MAX_SUBTYPE; ++i)
    if (fabs(pa.xi[0] - prev_pa.xi[i]) > eps)
      {
        cerr << "xi yet converged" << endl;
        return false;
      }

  for (int h=0; h<MAX_H; ++h)
    if (fabs(pa.omega[h] - prev_pa.omega[h]) > eps)
      {
        cerr << "omega yet converged" << endl;
        return false;
      }

  return true;
}

bool ancestor(int i, int q)
{
  while (i < q)
    q = calc_parent(q);

  if (i == q)
    return true;
  
  return false;
}


double calc_t(int i, std::vector<double> u)
{
  if (i == 1)
    return u[1];

  if (calc_remainder(i) == 0) // first child
    return u[i] * calc_t(calc_parent(i), u);

  return u[i] * calc_t(i - 1, u);
}

double calc_nu(int i, params* pa)
{
  // when MAX_CHILD = 3
  // if (i % MAX_CHILD == 2)
  //   return calc_nu(parent, beta) * (1 - beta[parent][0]) * beta[parent][1];

  // if (i % 3 == 0)
  //   return calc_nu(parent, beta) * (1 - beta[parent][0]) * (1 - beta[parent][1]) * beta[parent][2];

  // if (i % 3 == 1)
  //   return calc_nu(parent, beta) * (1 - beta[parent][0]) * (1 - beta[parent][1]) * (1 - beta[parent][2]);

  if (i == 0)
    return 1;

  if (i == 1)
    return 1 - pa->beta[0][0];
  
  int parent = calc_parent(i);
  int remainder = calc_remainder(i);
  double coef = 1;

  // cerr << "i = " << i << "\tparent = " << parent << "\tremainder = " << remainder << endl;
  
  for (int k=0; k<=remainder; ++k)
    coef *= 1 - pa->beta[parent][k];

  if (remainder < MAX_CHILD - 1)
    coef *= pa->beta[parent][remainder+1];
  
  return calc_nu(parent, pa) * coef;
}

double calc_n(int i, params* pa)
{
  double tmp = calc_nu(i, pa);
  
  // cerr << "calc_nu(" << i << ") = " << tmp << endl;
  if (i <= NONLEAF)
    // cerr << pa->beta[i][0] << endl;
  
  if (i <= NONLEAF)
    return tmp * (pa->beta)[i][0];
  return tmp;
}

void calc_n_all(std::vector<double>& n, params* pa)
{
  for (int i=0; i<=MAX_SUBTYPE; ++i)
    {
      double tmp = calc_n(i, pa);
      // cerr << "calc_n(" << i << ") = " << tmp << endl;
      n[i] = tmp;
    }
  
  return;
}

double calc_alpha(double t_i, double n_i)
{
  return log(CELL_MAX * n_i) / t_i;
}

double calc_mu(state st, std::vector<double> &n, double epsilon)
{
  // cout << "st.q = " << st.q << "\tst.h = " << st.h << endl;
  
  double numerator = 0;
  double denominator = 0;

  for (int i=0; i<=MAX_SUBTYPE; ++i)
    denominator += n[i] * st.l[i];

  for (int i=0; i<=MAX_SUBTYPE; ++i)
    {
      double x = 0;

      if (i == st.q)
        x = st.s;

      else
        {
          int prev_parent = i;
          int parent = calc_parent(i);
          while (parent > st.q)
            {
              prev_parent = parent;
              parent = calc_parent(prev_parent);
            }
          if (parent == st.q)
            {
              int remainder = (prev_parent + MAX_CHILD - 2) % MAX_CHILD;
              int h_j = st.h;
              for (int l = 0; l<remainder; ++l)
                h_j /= 2;
              h_j %= 2;
              if (h_j == 1)
                {
                  // cout << i << " is inherited" << endl;
                  x = 1;
                }
              
            }
        }
      
      numerator += n[i] * (epsilon * st.l[i] + (1 - 2 * epsilon) * x * st.r[i]);
    }
  
  return numerator / denominator;
}

double read_prob(READ re, state st, std::vector<double> &n, double epsilon)
{
  double mu = calc_mu(st, n, epsilon);

  // if (mu <= 0 || 1 <= mu)
  //   cout << "mu = " << mu << endl;
  
  // boost::math::binomial_distribution<double> binom(re.second, mu);

  double p = gsl_ran_binomial_pdf(re.first, mu, re.second);

  // if((boost::math::isnan)(p))
  //   {
  //     for (int i=0; i<=MAX_SUBTYPE; ++i)
  //       cout << st.l[i] << "\t" << st.r[i] << endl;
  //     cout << st.q << endl;
  //     cout << st.h << endl;
  //     cout << st.s << endl;
      
  //     cout << "pdf returned a NaN!" << endl;
  //     if (errno != 0)
  //       {
  //         // So errno has been set.
  //         cout << "errno is set to: " << errno << endl;
  //       }
  //   }
  // cout << "pdf" << tmp << endl;
  
  return p;
}

double variant_fraction(double x_i, int h, int i, double n_i, double t_i, double t_i_h, double beta_tilda)
{
  int s = (int) (x_i * FRACTIONS);

  return variant_fraction_normalized(s, h, i, n_i, t_i, t_i_h, beta_tilda);
  // return n_i / x_i / beta_tilda;
}

double d_t_variant_fraction(double x_i, int h, int i, double n_i, double t_i, double t_i_h, double beta_tilda)
{
  int s = (int) (x_i * FRACTIONS);
  double result, abserr;
  
  d_t_variant_fraction_normalized(s, h, i, n_i, t_i, t_i_h, beta_tilda, &result, &abserr);
  
  return result;
  // return 0;
}

double d_t_h_variant_fraction(double x_i, int h, int i, double n_i, double t_i, double t_i_h, double beta_tilda)
{
  int s = (int) (x_i * FRACTIONS);
  double result, abserr;
  
  d_t_h_variant_fraction_normalized(s, h, i, n_i, t_i, t_i_h, beta_tilda, &result, &abserr);
  
  return result;
  // return 0;
}

double d_n_variant_fraction(double x_i, int h, int i, double n_i, double t_i, double t_i_h, double beta_tilda)
{
  int s = (int) (x_i * FRACTIONS);
  double result, abserr;
  
  d_n_variant_fraction_normalized(s, h, i, n_i, t_i, t_i_h, beta_tilda, &result, &abserr);
  
  return result;
  // return 1.0 / x_i / beta_tilda;
}

int first_inherit_child(int i, int h)
{
  // if 0 is returned, no child inherit
  
  for (int j=0; j<MAX_CHILD; ++j)
    {
      if (h % 2 == 1)
        return calc_child(i, j);
      h /= 2;
    }
  return 0;
}

double llik(READS &res, states* sts, params* pa, hyperparams hpa)
{
  double sum = 0;
  
  // boost::math::beta_distribution<> beu(hpa.au, hpa.bu);
  // boost::math::beta_distribution<> bebeta(hpa.abeta, hpa.bbeta);
  
  for (int i=2; i<=MAX_SUBTYPE; ++i)
    sum += log(gsl_ran_beta_pdf(pa->u[i], hpa.au, hpa.bu));
    // sum += log(pdf(beu, pa->u[i]));

  sum += log(gsl_ran_beta_pdf(pa->beta[0][0], hpa.abeta, hpa.bbeta));
  // sum += log(pdf(bebeta, pa->beta[0][0]));
  
  for (int i=1; i <= NONLEAF; ++i)
    for (int j=0; j<MAX_CHILD; ++j)
      sum += log(gsl_ran_beta_pdf(pa->beta[i][j], hpa.abeta, hpa.bbeta));
      // sum += log(pdf(bebeta, pa->beta[i][j]));

  int count = 0;

  READS::iterator re = res.begin();
  
  for (READS::iterator re = res.begin(); re != res.end(); ++re)
    {
      for (states::iterator st = sts->begin(); st != sts->end(); ++st)
        {
          double a = 0;
      
          for (int i=1; i<=MAX_SUBTYPE; ++i)
            a += log(pa->pi[st->l[i]] * pa->kappa[st->l[i]][st->r[i]]);

          a += log(pa->xi[st->q] * pa->omega[st->h]);

          std::vector<double> n(MAX_SUBTYPE+1, 0);
          calc_n_all(n, pa);

          double t_q_h = -1;
          if (st->h != 0)
            {
              int ch = first_inherit_child(st->q, st->h);
              t_q_h = calc_t(ch, pa->u);
            }

          a += log(variant_fraction(st->s, st->h, st->q, n[st->q], calc_t(st->q, pa->u), t_q_h, pa->beta_tilda[st->q]) * read_prob(**re, *st, n, hpa.epsilon));

          sum += a;
        }
    }
  return sum;
}

typedef struct _rsph
{
  READS res;
  states* sts;
  V pi;
  VV kappa;
  V xi;
  V omega;
  V beta_tilda;
  hyperparams hpa;
}
  rsph;
  
double my_f(const gsl_vector *v, void *rsph_par)
{
  params* pa = new params;
  
  init_params(pa);

  for (int i=2; i<=MAX_SUBTYPE; ++i)
    pa->u[i] = gsl_vector_get(v, i-2);
  pa->beta[0][0] = gsl_vector_get(v, MAX_SUBTYPE-1);

  for (int i=1; i<=NONLEAF; ++i)
    for (int j=0; j<MAX_CHILD; ++j)
      pa->beta[i][j] = gsl_vector_get(v, MAX_SUBTYPE + MAX_CHILD*(i-1) + j);

  rsph* p = (rsph*) rsph_par;

  pa->pi = p->pi;
  pa->kappa = p->kappa;
  pa->xi = p->xi;
  pa->omega = p->omega;
  pa->beta_tilda = p->beta_tilda;
  
  return llik((p->res), p->sts, pa, p->hpa);
} // kokomade

double d_ln_betadist(double x, double a, double b)
{
  // zero if a = b
  
  // boost::math::beta_distribution<> beu(a, b);
  // boost::math::beta_distribution<> beu1(a - 1, b);
  // boost::math::beta_distribution<> beu2(a, b - 1);

  return (a + b - 1) * (gsl_ran_beta_pdf(x, a-1, b) - gsl_ran_beta_pdf(x, a, b-1)) / gsl_ran_beta_pdf(x, a, b);
  // return (a + b - 1) * (pdf(beu1, x) - pdf(beu2, x)) / pdf(beu, x);
}


double dllk_dui(states* sts, params* pa, hyperparams hpa, int i)
{
  double sum = 0;

  sum += d_ln_betadist(pa->u[i], hpa.au, hpa.bu);

  for (states::iterator st = sts->begin(); st != sts->end(); ++st)
    {
      int i_parent = calc_parent(i);
      int qk_parent = st->q;
      int qk_prev_parent = st->q;

      int q_hk = first_inherit_child(st->q, st->h);

      double tq = calc_t(st->q, pa->u);
      double tq_h = -1;
      if (st->h != 0)
        {
          tq_h = calc_t(q_hk, pa->u);
        }

      double nq = -1;
      double vf = -1;
      
      while (qk_parent > i_parent)
        {
          qk_prev_parent = qk_parent;
          qk_parent = calc_parent(qk_parent);
        }
      if (qk_parent == i_parent && i <= qk_prev_parent)
        {
          nq = calc_n(st->q, pa);
          vf =  variant_fraction(st->s, st->h, st->q, nq, tq, tq_h, pa->beta_tilda[st->q]);
          sum += tq * d_t_variant_fraction(st->s, st->h, st->q, nq, tq, tq_h, pa->beta_tilda[st->q]) / vf / pa->u[i];
        }

      if (st->h != 0)
        {
          int q_hk_parent = q_hk;
          int q_hk_prev_parent = q_hk;
      
          while (q_hk_parent > i_parent)
            {
              q_hk_prev_parent = q_hk_parent;
              q_hk_parent = calc_parent(q_hk_parent);
            }
          if (q_hk_parent == i_parent && i <= q_hk_prev_parent)
            {
              if (nq < 0)
                {
                  nq = calc_n(st->q, pa);
                  vf =  variant_fraction(st->s, st->h, st->q, nq, tq, tq_h, pa->beta_tilda[st->q]);
                }
              sum += tq_h * d_t_h_variant_fraction(st->s, st->h, st->q, nq, tq, tq_h, pa->beta_tilda[st->q]) / vf / pa->u[i];
            }
        }
    }
  return sum;
}

double dmu_dny(state st, std::vector<double> n, double epsilon, int y) // st is the state of y
{
  double numerator = 0;
  double denominator = 0;

  for (int i=0; i<MAX_SUBTYPE; ++i)
    denominator += n[i] * st.l[i];

  for (int i=0; i<MAX_SUBTYPE; ++i)
    {
      double x = 0;

      if (i == st.q)
        x = st.s;

      else
        {
          int prev_parent = i;
          int parent = calc_parent(i);
          while (parent > st.q)
            {
              prev_parent = parent;
              parent = calc_parent(prev_parent);
            }
          if (parent == st.q)
            {
              int remainder = calc_remainder(prev_parent);
              int h_j = st.h;
              for (int l = 0; l<remainder; ++l)
                h_j /= 2;
              h_j %= 2;
              if (h_j == 1)
                x = 1;
            }
        }
      numerator += n[i] * (epsilon * st.l[i] + (1 - 2*epsilon) * x * st.r[i]);
    }
  return (epsilon * st.l[y] + (1 - 2*epsilon) * st.s * st.r[y]) / denominator + st.l[y] * numerator / denominator / denominator;
}

double dreadprob_dmu(READ re, double mu)
{
  // boost::math::binomial_distribution<> binom(re.second - 1, mu);

  return re.second * (gsl_ran_binomial_pdf(re.first - 1, mu, re.second - 1) - gsl_ran_binomial_pdf(re.first, mu, re.second - 1));
  // return re.second * (pdf(binom, re.first - 1) - pdf(binom, re.first));
}

int child_order(int i, int q)
{
  int prev_parent = q;
  
  while (i < q)
    {
      prev_parent = q;
      q = calc_parent(q);
    }

  if (i == q)
    return calc_remainder(prev_parent);
  return -1;
}

double beta_factor(int y, int i, int j, params* pa)
{
  int l = -1;
  if ((l = child_order(i, y)) >= 0)
    {
      if (j == 0)
        {
          if (i == y)
            {
              if (i <= NONLEAF)
                return 1.0 / pa->beta[i][j];
            }
          else
            {
              return -1.0 / (1 - pa->beta[i][j]);
            }
        }
      else
        {
          if (j < l)
            return -1.0 / (1 - pa->beta[i][j]);
          else if (j == l)
            return 1.0 / pa->beta[i][j];
        }
    }
  return 0;
}

double dllik_dbeta(READS &res, states* sts, params* pa, hyperparams hpa, int i, int j)
{

  std::vector<double> n(MAX_SUBTYPE+1, 0);
  
  double sum = 0;

  sum += d_ln_betadist(pa->beta[i][j], hpa.abeta, hpa.bbeta);

  states::iterator st = sts->begin();
  for (READS::iterator re = res.begin(); re != res.end(); ++re)
    {
      double tq = calc_t(st->q, pa->u);
      int q_hk = first_inherit_child(st->q, st->h);
      double tq_h = -1;
      if (st->h != 0)
        tq_h = calc_t(q_hk, pa->u);

      calc_n_all(n, pa);

      double fac;
      
      sum += n[st->q] * beta_factor(st->q, i, j, pa) * d_n_variant_fraction(st->s, st->h, st->q, n[st->q], tq, tq_h, pa->beta_tilda[st->q]) / variant_fraction(st->s, st->h, st->q, n[st->q], tq, tq_h, pa->beta_tilda[st->q]);

      double acc = 0;

      double mu = calc_mu(*st, n, hpa.epsilon);
      
      for (int y=i; y<=MAX_SUBTYPE; ++y)
        {
          acc += n[y] * beta_factor(y, i, j, pa) * dmu_dny(*st, n, hpa.epsilon, y) * dreadprob_dmu(**re, mu);
        }
      sum += acc / read_prob(**re, *st, n, hpa.epsilon);
      
      ++st;
    }
  return sum;
}

void my_df(const gsl_vector *v, void *rsph_par, gsl_vector *df)
{
  params* pa = new params;
  
  init_params(pa);

  for (int i=2; i<=MAX_SUBTYPE; ++i)
    pa->u[i] = gsl_vector_get(v, i-2);
  pa->beta[0][0] = gsl_vector_get(v, MAX_SUBTYPE-1);

  for (int i=1; i<=NONLEAF; ++i)
    for (int j=0; j<MAX_CHILD; ++j)
      pa->beta[i][j] = gsl_vector_get(v, MAX_SUBTYPE + MAX_CHILD*(i-1) + j);

  rsph* p = (rsph*) rsph_par;

  pa->pi = p->pi;
  pa->kappa = p->kappa;
  pa->xi = p->xi;
  pa->omega = p->omega;
  pa->beta_tilda = p->beta_tilda;

  for (int i=2; i<=MAX_SUBTYPE; ++i)
    gsl_vector_set(df, i-2, dllk_dui(p->sts, pa, p->hpa, i));
  gsl_vector_set(df, MAX_SUBTYPE-1, dllik_dbeta((p->res), p->sts, pa, p->hpa, 0, 0));

  for (int i=1; i<=NONLEAF; ++i)
    for (int j=0; j<MAX_CHILD; ++j)
      gsl_vector_set(df, MAX_SUBTYPE + MAX_CHILD*(i-1) + j, dllik_dbeta((p->res), p->sts, pa, p->hpa, i, j));
}

void my_fdf(const gsl_vector *x, void *rsph_par, double *f, gsl_vector *df)
{
  *f = my_f(x, rsph_par);
  my_df(x, rsph_par, df);
}

double responsibility_numerator(state* st, params* pa_old, hyperparams hpa, READ re, std::vector<double> &n)
{
  double product = 1;
  
  for (int i=1; i<=MAX_SUBTYPE; ++i)
    product *= pa_old->pi[st->l[i]] * pa_old->kappa[st->l[i]][st->r[i]];
  
  int q_h = first_inherit_child(st->q, st->h);
  double tq_h = -1;
  if (st->h != 0)
    tq_h = calc_t(q_h, pa_old->u);

  product *= pa_old->xi[st->q] * pa_old->omega[st->h] * variant_fraction(st->s, st->h, st->q, n[st->q], calc_t(st->q, pa_old->u), tq_h, pa_old->beta_tilda[st->q]) * read_prob(re, *st, n, hpa.epsilon);
  
  return product;
}

double responsibility_partition_k(state* curr_st, params* pa_old, hyperparams hpa, READ re, std::vector<double> &n)
{
  double a = 0;
  
  for (curr_st->q = 1; curr_st->q <= MAX_SUBTYPE; ++curr_st->q)
    {
      // cout << "curr_st->q = " << curr_st->q << endl;
      
      double b = 0;

      for (curr_st->h = 0; curr_st->h < MAX_H; ++curr_st->h)
        {
          // cout << "curr_st->h = " << curr_st->h << endl;
          double c = 0;

          for (int s=1; s<=FRACTIONS; ++s)
            {
              curr_st->s = ((double)s) / FRACTIONS;
              // cout << "curr_st->s = " << curr_st->s << endl;

              int q_h = first_inherit_child(curr_st->q, curr_st->h);
              double tq_h = -1;
              if (curr_st->h != 0)
                {
                  tq_h = calc_t(q_h, pa_old->u);
                  // cout << "q_h: " << q_h << endl;
                  // cout << scientific;
                  // cout << "tq_h: " <<tq_h << endl;
                }

              c += variant_fraction(curr_st->s, curr_st->h, curr_st->q, n[curr_st->q], calc_t(curr_st->q, pa_old->u), tq_h, pa_old->beta_tilda[curr_st->q]) * read_prob(re, *curr_st, n, hpa.epsilon);
            }
          b += pa_old->omega[curr_st->h] * c;
        }
      a += pa_old->xi[curr_st->q] * b;
    }
  return a;
}

double responsibility_partition_subtype(state* curr_st, params* pa_old, hyperparams hpa, READ re, int curr_subtype, std::vector<double> &n)
{
  if (curr_subtype > MAX_SUBTYPE)
    return responsibility_partition_k(curr_st, pa_old, hpa, re, n);

  double t = 0;
  
  for (curr_st->l[curr_subtype] = 1; curr_st->l[curr_subtype] <= MAX_COPY; ++curr_st->l[curr_subtype])
    {
      // cout << "curr_st->l[" << curr_subtype << "] =" << curr_st->l[curr_subtype] << endl;
                                                         
      double s = 0;

      for (curr_st->r[curr_subtype] = 1; curr_st->r[curr_subtype] <= curr_st->l[curr_subtype]; ++curr_st->r[curr_subtype])
        {
          // cout << "curr_st->r[" << curr_subtype << "] =" << curr_st->r[curr_subtype] << endl;
          s += pa_old->kappa[curr_st->l[curr_subtype]][curr_st->r[curr_subtype]] * responsibility_partition_subtype(curr_st, pa_old, hpa, re, curr_subtype + 1, n);
        }
      // cout << endl;
      
      t += pa_old->pi[curr_st->l[curr_subtype]] * s;
    }
  return t;
}

double responsibility_partition(params* pa_old, hyperparams hpa, READ re, std::vector<double> &n)
{
  state* curr_st = new state;
  init_state(curr_st);

  curr_st->l[0] = 2;
  curr_st->r[0] = 0;
      
  return responsibility_partition_subtype(curr_st, pa_old, hpa, re, 1, n);
}

void responsibility_k(states* sts, state* curr_st, params* pa_old, hyperparams hpa, READ re, double resp_part, std::vector<double> &n)
{
  for (curr_st->q = 1; curr_st->q <= MAX_SUBTYPE; ++curr_st->q)
    {
      for (curr_st->h = 0; curr_st->h < MAX_H; ++curr_st->h)
        {
          for (int s=1; s<=FRACTIONS; ++s)
            {
              state* st = new state;
              init_state(st);

              st->k = curr_st->k;

              for (int i=1; i<=MAX_SUBTYPE; ++i)
                {
                  st->l[i] = curr_st->l[i];
                  st->r[i] = curr_st->r[i];
                }

              st->q = curr_st->q;
              st->h = curr_st->h;
              st->s = ((double)s) / FRACTIONS;

              st->resp = responsibility_numerator(st, pa_old, hpa, re, n) / resp_part;
              sts->push_back(*st);
            }
        }
    }
}

void responsibility_subtype(states* sts, state* curr_st, params* pa_old, hyperparams hpa, READ re, double resp_part, int curr_subtype, std::vector<double> &n)
{
  if (curr_subtype > MAX_SUBTYPE)
    {
      responsibility_k(sts, curr_st, pa_old, hpa, re, resp_part, n);
      return;
    }
  
    for (curr_st->l[curr_subtype] = 1; curr_st->l[curr_subtype] <= MAX_COPY; ++curr_st->l[curr_subtype])
    {
      for (curr_st->r[curr_subtype] = 1; curr_st->r[curr_subtype] <= curr_st->l[curr_subtype]; ++curr_st->r[curr_subtype])
        {
          responsibility_subtype(sts, curr_st, pa_old, hpa, re, resp_part, curr_subtype + 1, n);
        }
    }
  return;
}

void responsibility(states* sts, params* pa_old, hyperparams hpa, READS res)
{
  state* curr_st = new state;
  init_state(curr_st);
  
  curr_st->l[0] = 2;
  curr_st->r[0] = 0;

  std::vector<double> n(MAX_SUBTYPE+1, 0);
  calc_n_all(n, pa_old);

  for (READS::iterator re = res.begin(); re != res.end(); ++re)
    {
      double resp_part = responsibility_partition(pa_old, hpa, **re, n);
      
      responsibility_subtype(sts, curr_st, pa_old, hpa, **re, resp_part, 1, n);
    }
  return;
}

void maximization(params* pa_old, params* pa_new, hyperparams hpa, READS res, gsl_vector* x, int count, gsl_multimin_function_fdf my_func) // needs init x0
{
  states* sts = new states;
  responsibility(sts, pa_old, hpa, res);

  size_t iter = 0;
  int status;
  const gsl_multimin_fdfminimizer_type *T;
  gsl_multimin_fdfminimizer *s;
  
  if (count == 0)
    {
      x = gsl_vector_alloc(2*MAX_SUBTYPE - 1);
      for (int i=0; i<2*MAX_SUBTYPE-1; ++i)
        gsl_vector_set(x, i, 0.5);
    }

  rsph* p = new rsph;
  p->res = res;
  p->sts = sts;
  p->pi = pa_new->pi;
  p->kappa = pa_new->kappa;
  p->xi = pa_new->xi;
  p->omega = pa_new->omega;
  p->beta_tilda = pa_new->beta_tilda;
  p->hpa = hpa;

  my_func.params = p;

  T = gsl_multimin_fdfminimizer_conjugate_fr;
  s = gsl_multimin_fdfminimizer_alloc(T, 2*MAX_SUBTYPE-1);
  
  gsl_multimin_fdfminimizer_set(s, &my_func, x, 0.01, 1e-4);

  do
    {
      iter++;
      status = gsl_multimin_fdfminimizer_iterate(s);
      if (status)
        break;

      status = gsl_multimin_test_gradient(s->gradient, 1e-3);

      if (status == GSL_SUCCESS)
        printf ("Minimum found at:\n");

      cout << scientific;
      cout << "gsl_multimin_iter: " << (int)iter << endl;
      cout << "llik: " << s->f << endl;

      for (int i=2; i<=MAX_SUBTYPE; ++i)
        pa_new->u[i] = gsl_vector_get(s->x, i-2);
      pa_new->beta[0][0] = gsl_vector_get(s->x, MAX_SUBTYPE-1);
      
      for (int i=1; i<=NONLEAF; ++i)
        for (int j=0; j<MAX_CHILD; ++j)
          pa_new->beta[i][j] = gsl_vector_get(s->x, MAX_SUBTYPE + MAX_CHILD*(i-1) + j);
    } while (status == GSL_CONTINUE && iter < 10);
  gsl_multimin_fdfminimizer_free(s);
}



int main(int argc, char* argv[]) {
  _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK());
  gsl_set_error_handler_off ();
  
  cerr << MAX_SUBTYPE << endl;
  cerr << NONLEAF << endl;
  cerr << MAX_H << endl;
  
  ifstream infile;
  infile.open("reads.txt");

  READS* res = new READS;
  while(true)
    {
      READ* re = new READ;
      infile >> re->first >> re->second;
      if (infile.eof()) break;
      res->push_back(re);
    }
  infile.close();

  for (READS::iterator it=res->begin(); it!=res->end(); ++it)
    {
      cerr << (*it)->first << "\t" << (*it)->second << endl;
    }

  gegen.assign(FRACTIONS+1, V(GEGEN_MAX+1, 0));

  set_gegen(gegen);

  gegen_int.assign(FRACTIONS+1, 0);
  gegen_int_err.assign(FRACTIONS+1, 0);
  
  set_gegen_integral(gegen_int, gegen_int_err);
  
  gsl_multimin_function_fdf my_func;
  my_func.n = 2*MAX_SUBTYPE - 1;
  my_func.f = my_f;
  my_func.df = my_df;
  my_func.fdf = my_fdf;
  
  params* pa_old = new params; // pa_old initialization is needed
  init_params(pa_old);
  
  params* pa_new = new params;
  init_params(pa_new);
  
  hyperparams* hpa = new hyperparams; // hpa initialization is needed
  init_hyperparams(hpa);

  int em_maxit = 10;
  double em_eps = 1.0e-1;

  ofstream outfile;
  outfile.open(argv[1]);

  write_params(outfile, *pa_old);

  gsl_vector* x;
  for (int i=0; i<em_maxit; ++i)
    {
      cerr << "iter " << i << endl;
      outfile << "iter = " << i << endl;
      
      maximization(pa_old, pa_new, *hpa, *res, x, i, my_func); // res io is needed
      
      write_params(outfile, *pa_new);
      
      if (check_converge(*pa_new, *pa_old, em_eps)) // distance implementation is needed
        break;
      *pa_old = *pa_new;
      // delete pa_new;
    }
  delete_reads(res);
  outfile.close();

  return 0;
}
