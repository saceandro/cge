#include "em.hh"
#include "lbfgsb.hpp"
using namespace std;

bool ancestor(int i, int q)
{
  while (i < q)
    q = calc_parent(q);

  if (i == q)
    return true;
  
  return false;
}


double calc_t(int i, double* u)
{
  if (i == 1)
    return u[1];

  if ((i + J - 2) % J == 0) // first child
    return u[i] * calc_t(i/J + 1, u);

  return u[i] * calc_t(i - 1, u);
}

double calc_nu(int i, double beta[I][J])
{
  // when J = 3
  // if (i % J == 2)
  //   return calc_nu(parent, beta) * (1 - beta[parent][0]) * beta[parent][1];

  // if (i % 3 == 0)
  //   return calc_nu(parent, beta) * (1 - beta[parent][0]) * (1 - beta[parent][1]) * beta[parent][2];

  // if (i % 3 == 1)
  //   return calc_nu(parent, beta) * (1 - beta[parent][0]) * (1 - beta[parent][1]) * (1 - beta[parent][2]);

  if (i == 0)
    return 1;

  int parent = (i + J - 2) / J;
  int remainder = parent % J;
  double coef = 1;

  for (int k=0; k<=remainder; ++k)
    coef *= 1 - beta[parent][k];

  if (remainder < J - 1)
    coef *= beta[parent][remainder+1];
  
  return calc_nu(parent, beta) * coef;
}

double calc_n(int i, double* beta[J])
{
  return calc_nu(i, beta) * beta[i][0];
}

void calc_n_all(double* n, double* beta[J])
{
  for (int i=0; i<=I; ++i)
    n[i] = calc_n(i, beta);
  return;
}

double calc_alpha(double t_i, double n_i, int N)
{
  return log(N * n_i) / t_i;
}

double calc_mu(state st, double* n, double epsilon)
{
  double numerator = 0;
  double denominator = 0;

  for (int i=0; i<I; ++i)
    denominator += n[i] * st.l[i];

  for (int i=0; i<I; ++i)
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
              int remainder = (prev_parent + J - 2) % J;
              int h_j = st.h;
              for (int l = 0; l<remainder; ++l)
                h_j /= 2;
              h_j %= 2;
              if (h_j == 1)
                x = 1;
            }
        }
      
      numerator += n[i] * (epsilon * st.l[i] + (1 - 2 * epsilon) * x * st.r[i]);
    }
  
  double mu = numerator / denominator;
}

double read_prob(read re, state st, double* n, double epsilon)
{
  double mu = calc_mu(st, n, epsilon);

  boost::math::binomial_distribution<> binom(re.first, mu);

  return pdf(binom, re.second);
}



double variant_fraction(double x_i, int h, int i, double n_i, double t_i, double beta_tilda)
{
  return n_i / x_i / beta_tilda;
}

double d_t_variant_fraction(double x_i, int h, int i, double n_i, double t_i, double beta_tilda)
{
  return 0;
}

double d_n_variant_fraction(double x_i, int h, int i, double n_i, double t_i, double beta_tilda)
{
  return 1.0 / x_i / beta_tilda;
}

// double variant_fraction_alpha(double x_i, int h, int i, double alpha_i, double t_i, double beta_tilda, double N)
// {
//   return exp(alpha_i * t_i)  / x_i / beta_tilda / N;
// }

// double d_alpha_variant_fraction(double x_i, int h, int i, double alpha_i, double t_i, double beta)
// {
//   return t_i * variant_fraction(x_i, h, i, alpha_i, t_i, beta);
// }

// double d_t_variant_fraction(double x_i, int h, int i, double alpha_i, double t_i, double beta)
// {
//   return alpha_i * variant_fraction(x_i, h, i, alpha_i, t_i, beta);
// }

double llik(reads res, states sts, params pa, hyperparams hpa)
{
  double sum = 0;
  
  boost::math::beta_distribution<> beu(hpa.au, hpa.bu);
  boost::math::beta_distribution<> bebeta(hpa.abeta, hpa.bbeta);
  
  for (int i=2; i<=I; ++i)
    sum += log(pdf(beu, params.u[i]));

  sum += log(pdf(bebeta, params.beta[0][0]));
  
  for (int i=1; i <= I-J**(D-1); ++i)
    for (int j=0; j<J; ++j)
      sum += log(pdf(bebeta, params.beta[i][j]));

  for (state st = sts.begin(), std::vector<read>::iterator re = res.begin(); re != res.end(); ++sts, ++re)
    {
      double a = 0;
      
      for (int i=1; i<=I; ++i)
        a += log(params.pi[st.l[i]] * params.kappa[st.l[i]][st.r[i]]);

      a += log(params.xi[st.q] * params.omega[st.h]);

      double n[I+1];
      calc_n_all(n, pa.beta);
      a += log(variant_fraction(st.s, st.h, st.q, n, calc_t(st.q, pa.u), pa.beta_tilda[st.q]) * read_prob(re, st, n, hpa.epsilon));

      sum += a;
    }
  return sum;
}

double d_ln_betadist(double x, double a, double b)
{
  boost::math::beta_distribution<> beu(a, b);
  boost::math::beta_distribution<> beu1(a - 1, b);
  boost::math::beta_distribution<> beu2(a, b - 1);
  
  return (a + b - 1) * (pdf(beu1, x) - pdf(beu2, x)) / pdf(beu, x);
}


double dllk_dui(states sts, params pa, hyperparams hpa, int i)
{
  double sum = 0;

  sum += d_ln_betadist(pa.u[i], hpa.au, hpa.bu);

  for (state st = sts.begin(); st != sts.end(); ++sts)
    {
      int i_parent = calc_parent(i);
      int qk_paarent = st.q;
      int qk_prev_parent = st.q;
      
      while (qk_parent > i_parent)
        {
          qk_prev_parent = qk_parent;
          qk_parent = calc_parent(qk_parent);
        }
      if (qk_parent == i_parent && i <= qk_prev_parent)
        {
          double tq = calc_t(st.q, pa.u);
          double nq = calc_n(st.q, pa.beta);
          sum += tq * d_t_variant_fraction(st.s, st.h, st.q, nq, nt, pa.beta_tilda[st.q]) / variant_fraction(st.s, st.h, st.q, nq, tq, pa.beta_tilda[st.q]) / pa.u[i];
        }
    }
  return sum;
}

double dmu_dny(state st, double* n, double epsilon, int y) // st is the state of y
{
  double numerator = 0;
  double denominator = 0;

  for (int i=0; i<I; ++i)
    denominator += n[i] * st.l[i];

  for (int i=0; i<I; ++i)
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
  return (epsilon * st.l[y] + (1 - 2*epsilon) * st.s * st.r[y]) / denominator + st.l[y] * numerator / denominator**2;
}

double dreadprob_dmu(read re, double mu)
{
  boost::math::binomial_distribution<> binom(re.first - 1, mu);

  return re.first * (pdf(binom, re.second - 1) - pdf(binom, re.second));
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

doubke dllik_dbeta(reads res, states sts, params pa, hyperparams hpa, int i, int j)
{

  double n[I+1];
  
  double sum = 0;

  sum += d_ln_betadist(pa.beta[i][j], hpa.abeta, hpa.bbeta);

  for (read re = res.begin(), state st = sts.begin(); re != res.end(); ++re, ++sts)
    {
      double tq = calc_t(st.q, pa.u);

      calc_n_all(n, pa.beta);

      sum += d_n_variant_fraction(st.s, st.h, st.q, n[st.q], nt, pa.beta_tilda[st.q]) / variant_fraction(st.s, st.h, st.q, n[st.q], tq, pa.beta_tilda[st.q]);

      double acc = 0;

      calc_mu(st, n, epsilon);
      
      for (int y=i; y<=I; ++y)
        {
          double factor = n[y] * dmu_dny(st, n, hpa.epsilon, y) * dreadprob_dmu(re, mu);

          int l = -1;
          if ((l = child_order(i, y)) >= 0)
            {
              if (j == 0)
                {
                  if (i == y)
                    {
                      if (i <= NONLEAF)
                        acc += factor / pa.beta[i][j];
                    }
                  else
                    {
                      acc -= factor / (1 - pa.beta[i][j])
                    }
                }
              else
                {
                  if (j < l)
                    acc -= factor / (1 - pa.beta[i][j]);
                  else if (j == l)
                    acc += factor / pa.beta[i][j];
                }
            }

        }
      sum += acc / read_prob(re, st, n, hpa.epsilon);
    }
  return sum;
}


double responsibility_numerator(state st, params pa_old, hyperparams hpa, read re, double* n)
{
  double product = 1;
  
  for (int i=1; i<I; i++)
    product *= pa_old.pi[st.l[i]] * pa_old.kappa[st.l[i]][st.r[i]];
  product *= pa_old.xi[st.q] * pa_old.omega[st.h] * variant_fraction(st.s, st.h, st.q, n[st.q], calc_t(st.q, pa_old.u), pa_old.beta_tilda[st.q]) * read_prob(re, st, n, hpa.epsilon);
  
  return product;
}

double responsibility_partition_k(state curr_st, params pa_old, hyperparams hpa, read re, double* n)
{
  double a = 0;
  
  for (curr_st.q = 1; curr_st.q < I; ++curr_st.q)
    {
      double b = 0;

      for (curr_st.h = 0; curr_st.h < MAX_H; ++curr_st.h)
        {
          double c = 0;

          for (int s=0; s<FRACTIONS; ++s)
            {
              curr_st.s = ((double)s) / FRACTIONS;
              c += variant_fraction(curr_st.s, curr_st.h, curr_st.q, n[curr_st.q], calc_t(curr_st.q, pa_old.u), pa_old.beta_tilda[curr_st.q]) * read_prob(re, curr_st, hpa.epsilon);
            }
          b += pa_old.omega[curr_st.h] * c;
        }
      a += pa_old.xi[curr_st.q] * b;
    }
  return a;
}

double responsibility_partition_subtype(state curr_st, params pa_old, hyperparams hpa, read re, int curr_subtype, double* n)
{
  if (curr_subtype > I)
    return responsibility_partition_k(curr_st, pa_old, hpa, re, n);

  double t = 0;
  
  for (curr_st.l[curr_subtype] = 0; curr_st.l[curr_subtype] < C; ++curr_st.l[curr_subtype])
    {
      double s = 0;

      for (curr_st.r[curr_subtype] = 0; curr_st.r[curr_subtype] < curr_st.l[curr_subtype]; ++curr_st.r[curr_subtype])
        {
          s += pa_old.kappa[curr_st.l[curr_subtype]][curr_st.r[curr_subtype]] * responsibility_partition_subtype(curr_st, pa_old, hpa, re, ++curr_subtype, k, n);
        }
      t += pa_old.pi[curr_st.l[curr_subtype]] * s;
    }
  return t;
}

double responsibility_partition(params pa_old, hyperparams hpa, read re, double* n)
{
  state curr_st;
  
  curr_st.l[0] = 2;
  curr_st.r[0] = 0;
      
  return responsilility_partition_subtype(curr_st, pa_old, hpa, re, 1, n);
}

void responsibility_k(states sts, state curr_st, params pa_old, hyperparams hpa, read re, double resp_part, double* n)
{
  for (curr_st.q = 1; curr_st.q < I; ++curr_st.q)
    {
      for (curr_st.h = 0; curr_st.h < MAX_H; ++curr_st.h)
        {
          for (int s=0; s<FRACTIONS; ++s)
            {
              state st = new state;
              
              st.k = curr_st.k;

              for (int i=1; i<=I; ++i)
                {
                  st.l[i] = curr_st.l[i];
                  st.r[i] = curr_st.r[i];
                }

              st.q = curr_st.q;
              st.h = curr_st.h;
              st.s = ((double)s) / FRACTIONS;

              st.resp = responsibility_numerator(st, pa_old, hpa, re, n) / resp_part;
              sts.push_back(st);
            }
        }
    }
}

void responsilibity_subtype(states sts, state curr_st, params pa_old, hyperparams hpa, read re, double resp_part, double* n)
{
  if (curr_subtype > I)
    return;

  for (curr_st.l[curr_subtype] = 0; curr_st.l[curr_subtype] < C; ++curr_st.l[curr_subtype])
    {
      for (curr_st.r[curr_subtype] = 0; curr_st.r[curr_subtype] < curr_st.l[curr_subtype]; ++curr_st.r[curr_subtype])
        {
          responsibility_k(sts, curr_st, pa_old, hpa, re, resp_part, n);
        }
    }
  return;
}

void responsibility(states sts, params pa_old, hyperparams hpa, reads res)
{
  state curr_st;

  curr_st.l[0] = 2;
  curr_st.r[0] = 0;

  double* n[I+1];
  calc_n_all(n, pa_old.beta);

  for (curr_st.k=0, std::vector<read>::iterator re = res.begin(); re != res.end(); ++k, ++re)
    {
      double resp_part = responsibility_partition(old_pa, re, n);
      responsibility_subtype(sts, curr_st, pa_old, hpa, re, resp_part, n);
    }
  return;
}

struct ComputePdubeta {
  ComputePdubeta(reads res, states sts, params pa, hyperparams hpa) : _res(res), _sts(sts), _pa(pa), _hpa(hpa), _count(0) {}
  int operator()(const double& x, double& fn, std::vector<double>& gr) {
    for (int i=2; i<=I; ++i)
      _pa.u[i] = x[i-2];
    _pa.beta[0][0] = x[I-1];

    for (int i=1; i<=NONREAF; ++i)
      for (int j=0; j<J; ++j)
        _pa.beta[i][j] = x[I + J*(i-1) + j];
    
    fn = llik(_res, _sts, _pa, _hpa)

    gr.assign(x.size(), 0);
    
    for (int i=2; i<=I; ++i)
      gr[i-2] = dllk_dui(_sts, _pa, _hpa, i);
    gr[I-1] = dllik_dbeta(_res, _sts, _pa, _hpa, 0, 0);
    
    for (int i=1; i<=NONREAF; ++i)
      for (int j=0; j<J; ++j)
        gr[I + J*(i-1) + j] = dllik_dbeta(_res, _sts, _pa, _hpa, i, j);
    
    cout << "i:" << ++_count << ",fn:" << fn << endl;
    return 0;
  }
  
  reads _res;
  states _sts;
  params _pa;
  hyperparams _hpa;
  int    _count;
};

void maximization(params pa_old, params new_pa, hyperparams hya, reads res, Lbfgsb minimizer)
{
  states sts = new states;
  responsibility(sts, pa_old, hpa, res);

  double partition;
  for (std::vector<state>::iterator st = sts.begin(); st != sts.end(); ++st)
    partition += st.resp;

  double pi_numerator[C+1];
  for (int l=1; l<=C; ++l)
    {
      for (std::vector<state>::iterator st = sts.begin(); st != sts.end(); ++st)
        {
          int count = 0;
          for (int i=1; i<=I; ++i)
            {
              if (st.l[i] == l)
                count++;
            }
          pi_numerator[l] += count * st.resp;
        }
      new_pa.pi[l] = pi_numerator[l] / partition / I;
    }

  double kappa_numerator[C][C];
  for (int l=1; l<=C; ++l)
    {
      for (int r=1; r<=l; ++r)
        {
          for (std::vector<state>::iterator st = sts.begin(); st != sts.end(); ++st)
            {
              int count = 0;
              for (int i=1; i<=I; ++i)
                {
                  if (st.l[i] == l && st.r[i] == r)
                    count++;
                }
              kappa_numerator[l][r] += count * st.resp;
            }
          new_pa.kappa[l][r] = kappa_numerator[l][r] / pi_numerator[l];
        }
    }

  double xi_numerator[I+1];
  for (int q=1; q<=I; ++q)
    {
      for (std::vector<state>::iterator st = sts.begin(); st != sts.end(); ++st)
        {
          if (st.q == q)
            xi_numerator[q] += st.resp;
        }
      new_pa.xi[q] = xi_numerator[q] / partition;
    }

  double omega_numerator[MAX_H];
  for (int h=0; h<MAX_H; ++q)
    {
      for (std::vector<state>::iterator st = sts.begin(); st != sts.end(); ++st)
        {
          if (st.h == h)
            xi_numerator[h] += st.resp;
        }
      new_pa.omega[h] = omega_numerator[h] / partition;
    }

  minimizer.minimize(minimizer.best_x(), ComputePdubeta(res, sts, pa, hpa));
  cout << "iter:" << minimizer.iter() << ",best_fn:" << minimizer.best_fn() << ",best_x:";
  const vector<double>& x = minimizer.best_x();
  for (int i = 0; i < (int)x.size(); ++i) { cout << x[i] << ",";}
  cout << endl;

  for (int i=2; i<=I; ++i)
    new_pa.u[i] = x[i-2];
  new_pa.beta[0][0] = x[I-1];
  
  for (int i=1; i<=NONREAF; ++i)
    for (int j=0; j<J; ++j)
      new_pa.beta[i][j] = x[I + J*(i-1) + j];
}

int main(int argc, char* argv[]) {
  double c = 1;
  Lbfgsb minimizer;
  minimizer.set_eps(1.0e-4);
  minimizer.set_maxit(30);
  vector<double> x0(2I-1, 0.5);

  int em_maxit = 10;
  int em_eps = 1.0e-1;
  
  int iter = 0;
  params pa_old = new params; // pa_old initialization is needed
  params pa_new = new params;
  hyperparams hpa = new hyperparams; // hpa initialization is needed
  
  for (int i=0; i<em_maxit; ++i)
    {
      maximization(pa_old, new_pa, hpa, res, minimizer); // res io is needed
      if (distance(pa_old, pa_new) < em_eps) // distance implementation is needed
        break;
      delete pa_old;
      pa_old = pa_new;
    }
  return 0;
}
