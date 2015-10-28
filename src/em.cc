#include "em.hh"
#include "lbfgsb.hpp"
using namespace std;

// params::params(std::vector<double> _u, std::vector<std::vector<double> > _beta, double* _pi, double** _kappa, std::vector<double> _xi, std::vector<double> _omega, std::vector<double> _beta_tilda)
// {
// }

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

  if ((i + MAX_CHILD - 2) % MAX_CHILD == 0) // first child
    return u[i] * calc_t(i/MAX_CHILD + 1, u);

  return u[i] * calc_t(i - 1, u);
}

double calc_nu(int i, std::vector<std::vector<double> > beta)
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

  int parent = (i + MAX_CHILD - 2) / MAX_CHILD;
  int remainder = parent % MAX_CHILD;
  double coef = 1;

  for (int k=0; k<=remainder; ++k)
    coef *= 1 - beta[parent][k];

  if (remainder < MAX_CHILD - 1)
    coef *= beta[parent][remainder+1];
  
  return calc_nu(parent, beta) * coef;
}

double calc_n(int i, std::vector<std::vector<double> > beta)
{
  return calc_nu(i, beta) * beta[i][0];
}

void calc_n_all(std::vector<double> n, std::vector<std::vector<double> > beta)
{
  for (int i=0; i<=MAX_SUBTYPE; ++i)
    n[i] = calc_n(i, beta);
  return;
}

double calc_alpha(double t_i, double n_i, int N)
{
  return log(N * n_i) / t_i;
}

double calc_mu(state st, std::vector<double> n, double epsilon)
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
              int remainder = (prev_parent + MAX_CHILD - 2) % MAX_CHILD;
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
  
  return numerator / denominator;
}

double read_prob(READ re, state st, std::vector<double> n, double epsilon)
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

double llik(READS res, states sts, params pa, hyperparams hpa)
{
  double sum = 0;
  
  boost::math::beta_distribution<> beu(hpa.au, hpa.bu);
  boost::math::beta_distribution<> bebeta(hpa.abeta, hpa.bbeta);
  
  for (int i=2; i<=MAX_SUBTYPE; ++i)
    sum += log(pdf(beu, pa.u[i]));

  sum += log(pdf(bebeta, pa.beta[0][0]));
  
  for (int i=1; i <= NONLEAF; ++i)
    for (int j=0; j<MAX_CHILD; ++j)
      sum += log(pdf(bebeta, pa.beta[i][j]));

  int count = 0;
  for (std::vector<state>::iterator st = sts.begin(); st != sts.end(); ++st)
    {
      READ re = res[count];
      double a = 0;
      
      for (int i=1; i<=MAX_SUBTYPE; ++i)
        a += log(pa.pi[st->l[i]] * pa.kappa[st->l[i]][st->r[i]]);

      a += log(pa.xi[st->q] * pa.omega[st->h]);

      std::vector<double> n(MAX_SUBTYPE+1, 0);
      calc_n_all(n, pa.beta);
      a += log(variant_fraction(st->s, st->h, st->q, n[st->q], calc_t(st->q, pa.u), pa.beta_tilda[st->q]) * read_prob(re, *st, n, hpa.epsilon));

      sum += a;
      count++;
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

  for (std::vector<state>::iterator st = sts.begin(); st != sts.end(); ++st)
    {
      int i_parent = calc_parent(i);
      int qk_parent = st->q;
      int qk_prev_parent = st->q;
      
      while (qk_parent > i_parent)
        {
          qk_prev_parent = qk_parent;
          qk_parent = calc_parent(qk_parent);
        }
      if (qk_parent == i_parent && i <= qk_prev_parent)
        {
          double tq = calc_t(st->q, pa.u);
          double nq = calc_n(st->q, pa.beta);
          sum += tq * d_t_variant_fraction(st->s, st->h, st->q, nq, tq, pa.beta_tilda[st->q]) / variant_fraction(st->s, st->h, st->q, nq, tq, pa.beta_tilda[st->q]) / pa.u[i];
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

double dllik_dbeta(READS res, states sts, params pa, hyperparams hpa, int i, int j)
{

  std::vector<double> n(MAX_SUBTYPE+1, 0);
  
  double sum = 0;

  sum += d_ln_betadist(pa.beta[i][j], hpa.abeta, hpa.bbeta);

  states::iterator st = sts.begin();
  for (READS::iterator re = res.begin(); re != res.end(); ++re)
    {
      double tq = calc_t(st->q, pa.u);

      calc_n_all(n, pa.beta);

      sum += d_n_variant_fraction(st->s, st->h, st->q, n[st->q], tq, pa.beta_tilda[st->q]) / variant_fraction(st->s, st->h, st->q, n[st->q], tq, pa.beta_tilda[st->q]);

      double acc = 0;

      double mu = calc_mu(*st, n, hpa.epsilon);
      
      for (int y=i; y<=MAX_SUBTYPE; ++y)
        {
          double factor = n[y] * dmu_dny(*st, n, hpa.epsilon, y) * dreadprob_dmu(*re, mu);

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
                      acc -= factor / (1 - pa.beta[i][j]);
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
      sum += acc / read_prob(*re, *st, n, hpa.epsilon);
      ++st;
    }
  return sum;
}


double responsibility_numerator(state st, params pa_old, hyperparams hpa, READ re, std::vector<double> n)
{
  double product = 1;
  
  for (int i=1; i<MAX_SUBTYPE; i++)
    product *= pa_old.pi[st.l[i]] * pa_old.kappa[st.l[i]][st.r[i]];
  product *= pa_old.xi[st.q] * pa_old.omega[st.h] * variant_fraction(st.s, st.h, st.q, n[st.q], calc_t(st.q, pa_old.u), pa_old.beta_tilda[st.q]) * read_prob(re, st, n, hpa.epsilon);
  
  return product;
}

double responsibility_partition_k(state curr_st, params pa_old, hyperparams hpa, READ re, std::vector<double> n)
{
  double a = 0;
  
  for (curr_st.q = 1; curr_st.q < MAX_SUBTYPE; ++curr_st.q)
    {
      double b = 0;

      for (curr_st.h = 0; curr_st.h < MAX_H; ++curr_st.h)
        {
          double c = 0;

          for (int s=0; s<FRACTIONS; ++s)
            {
              curr_st.s = ((double)s) / FRACTIONS;
              c += variant_fraction(curr_st.s, curr_st.h, curr_st.q, n[curr_st.q], calc_t(curr_st.q, pa_old.u), pa_old.beta_tilda[curr_st.q]) * read_prob(re, curr_st, n, hpa.epsilon);
            }
          b += pa_old.omega[curr_st.h] * c;
        }
      a += pa_old.xi[curr_st.q] * b;
    }
  return a;
}

double responsibility_partition_subtype(state curr_st, params pa_old, hyperparams hpa, READ re, int curr_subtype, std::vector<double> n)
{
  if (curr_subtype > MAX_SUBTYPE)
    return responsibility_partition_k(curr_st, pa_old, hpa, re, n);

  double t = 0;
  
  for (curr_st.l[curr_subtype] = 0; curr_st.l[curr_subtype] < MAX_COPY; ++curr_st.l[curr_subtype])
    {
      double s = 0;

      for (curr_st.r[curr_subtype] = 0; curr_st.r[curr_subtype] < curr_st.l[curr_subtype]; ++curr_st.r[curr_subtype])
        {
          s += pa_old.kappa[curr_st.l[curr_subtype]][curr_st.r[curr_subtype]] * responsibility_partition_subtype(curr_st, pa_old, hpa, re, curr_subtype, n);
          curr_subtype++;
        }
      t += pa_old.pi[curr_st.l[curr_subtype]] * s;
    }
  return t;
}

double responsibility_partition(params pa_old, hyperparams hpa, READ re, std::vector<double> n)
{
  state curr_st;
  
  curr_st.l[0] = 2;
  curr_st.r[0] = 0;
      
  return responsibility_partition_subtype(curr_st, pa_old, hpa, re, 1, n);
}

void responsibility_k(states sts, state curr_st, params pa_old, hyperparams hpa, READ re, double resp_part, std::vector<double> n)
{
  for (curr_st.q = 1; curr_st.q < MAX_SUBTYPE; ++curr_st.q)
    {
      for (curr_st.h = 0; curr_st.h < MAX_H; ++curr_st.h)
        {
          for (int s=0; s<FRACTIONS; ++s)
            {
              state* st = new state;
              
              st->k = curr_st.k;

              for (int i=1; i<=MAX_SUBTYPE; ++i)
                {
                  st->l[i] = curr_st.l[i];
                  st->r[i] = curr_st.r[i];
                }

              st->q = curr_st.q;
              st->h = curr_st.h;
              st->s = ((double)s) / FRACTIONS;

              st->resp = responsibility_numerator(*st, pa_old, hpa, re, n) / resp_part;
              sts.push_back(*st);
            }
        }
    }
}

void responsibility_subtype(states sts, state curr_st, params pa_old, hyperparams hpa, READ re, double resp_part, int curr_subtype, std::vector<double> n)
{
  if (curr_subtype > MAX_SUBTYPE)
    return;

  for (curr_st.l[curr_subtype] = 0; curr_st.l[curr_subtype] < MAX_COPY; ++curr_st.l[curr_subtype])
    {
      for (curr_st.r[curr_subtype] = 0; curr_st.r[curr_subtype] < curr_st.l[curr_subtype]; ++curr_st.r[curr_subtype])
        {
          responsibility_subtype(sts, curr_st, pa_old, hpa, re, resp_part, curr_subtype++, n);
        }
    }
  return;
}

void responsibility(states sts, params pa_old, hyperparams hpa, READS res)
{
  state curr_st;

  curr_st.l[0] = 2;
  curr_st.r[0] = 0;

  std::vector<double> n(MAX_SUBTYPE+1, 0);
  calc_n_all(n, pa_old.beta);

  for (READS::iterator re = res.begin(); re != res.end(); ++re)
    {
      double resp_part = responsibility_partition(pa_old, hpa, *re, n);
      
      responsibility_subtype(sts, curr_st, pa_old, hpa, *re, resp_part, 1, n);
    }
  return;
}

struct ComputePdubeta {
  ComputePdubeta(READS res, states sts, params pa, hyperparams hpa) : _res(res), _sts(sts), _pa(pa), _hpa(hpa), _count(0) {}
  int operator()(const std::vector<double> x, double& fn, std::vector<double>& gr) { // x: reference?
    for (int i=2; i<=MAX_SUBTYPE; ++i)
      _pa.u[i] = x[i-2];
    _pa.beta[0][0] = x[MAX_SUBTYPE-1];

    for (int i=1; i<=NONLEAF; ++i)
      for (int j=0; j<MAX_CHILD; ++j)
        _pa.beta[i][j] = x[MAX_SUBTYPE + MAX_CHILD*(i-1) + j];
    
    fn = llik(_res, _sts, _pa, _hpa);

    gr.assign(x.size(), 0);
    
    for (int i=2; i<=MAX_SUBTYPE; ++i)
      gr[i-2] = dllk_dui(_sts, _pa, _hpa, i);
    gr[MAX_SUBTYPE-1] = dllik_dbeta(_res, _sts, _pa, _hpa, 0, 0);
    
    for (int i=1; i<=NONLEAF; ++i)
      for (int j=0; j<MAX_CHILD; ++j)
        gr[MAX_SUBTYPE + MAX_CHILD*(i-1) + j] = dllik_dbeta(_res, _sts, _pa, _hpa, i, j);
    
    cout << "i:" << ++_count << ",fn:" << fn << endl;
    return 0;
  }
  
  READS _res;
  states _sts;
  params _pa;
  hyperparams _hpa;
  int    _count;
};

void maximization(params pa_old, params new_pa, hyperparams hpa, READS res, Lbfgsb minimizer)
{
  states sts;
  responsibility(sts, pa_old, hpa, res);

  double partition = 0;
  for (std::vector<state>::iterator st = sts.begin(); st != sts.end(); ++st)
    partition += st->resp;

  double pi_numerator[MAX_COPY+1];
  for (int l=1; l<=MAX_COPY; ++l)
    {
      for (std::vector<state>::iterator st = sts.begin(); st != sts.end(); ++st)
        {
          int count = 0;
          for (int i=1; i<=MAX_SUBTYPE; ++i)
            {
              if (st->l[i] == l)
                count++;
            }
          pi_numerator[l] += count * st->resp;
        }
      new_pa.pi[l] = pi_numerator[l] / partition / MAX_SUBTYPE;
    }

  double kappa_numerator[MAX_COPY][MAX_COPY];
  for (int l=1; l<=MAX_COPY; ++l)
    {
      for (int r=1; r<=l; ++r)
        {
          for (std::vector<state>::iterator st = sts.begin(); st != sts.end(); ++st)
            {
              int count = 0;
              for (int i=1; i<=MAX_SUBTYPE; ++i)
                {
                  if (st->l[i] == l && st->r[i] == r)
                    count++;
                }
              kappa_numerator[l][r] += count * st->resp;
            }
          new_pa.kappa[l][r] = kappa_numerator[l][r] / pi_numerator[l];
        }
    }

  double xi_numerator[MAX_SUBTYPE+1];
  for (int q=1; q<=MAX_SUBTYPE; ++q)
    {
      for (std::vector<state>::iterator st = sts.begin(); st != sts.end(); ++st)
        {
          if (st->q == q)
            xi_numerator[q] += st->resp;
        }
      new_pa.xi[q] = xi_numerator[q] / partition;
    }

  double omega_numerator[MAX_H];
  for (int h=0; h<MAX_H; ++h)
    {
      for (std::vector<state>::iterator st = sts.begin(); st != sts.end(); ++st)
        {
          if (st->h == h)
            xi_numerator[h] += st->resp;
        }
      new_pa.omega[h] = omega_numerator[h] / partition;
    }

  minimizer.minimize(minimizer.best_x(), ComputePdubeta(res, sts, new_pa, hpa));
  cout << "iter:" << minimizer.iter() << ",best_fn:" << minimizer.best_fn() << ",best_x:";
  const vector<double>& x = minimizer.best_x();
  for (int i = 0; i < (int)x.size(); ++i) { cout << x[i] << ",";}
  cout << endl;

  for (int i=2; i<=MAX_SUBTYPE; ++i)
    new_pa.u[i] = x[i-2];
  new_pa.beta[0][0] = x[MAX_SUBTYPE-1];
  
  for (int i=1; i<=NONLEAF; ++i)
    for (int j=0; j<MAX_CHILD; ++j)
      new_pa.beta[i][j] = x[MAX_SUBTYPE + MAX_CHILD*(i-1) + j];
}

void init_params(params pa)
{
  for (int i=2; i<=MAX_SUBTYPE; ++i)
    pa.u[i] = 0.5;

  pa.beta[0][0] = 0.5;
  for (int i=1; i<=NONLEAF; ++i)
    for (int j=0; j<MAX_CHILD; ++j)
      pa.beta[i][j] = 0.5;

  for (int c=1; c<=MAX_COPY; ++c)
    pa.pi[c] = 1.0 / MAX_COPY;

  for (int c=1; c<=MAX_COPY; ++c)
    for (int d=1; d<=c; ++d)
      pa.kappa[c][d] = 1.0 / c;

  for (int i=1; i<=MAX_SUBTYPE; ++i)
    pa.xi[i] = 1.0 / MAX_SUBTYPE;

  for (int h=0; h<MAX_H; ++h)
    pa.omega[h] = 1.0 / MAX_H;

  for (int i=1; i<=MAX_SUBTYPE; ++i)
    pa.beta_tilda[i] = 0.5;

  return;
}

void init_hyperparams(hyperparams hpa)
{
  hpa.au = 2;
  hpa.bu = 2;
  hpa.abeta = 2;
  hpa.bbeta = 2;
  hpa.epsilon = 0.05;
  
  return;
}

bool check_converge(params pa, params prev_pa,  double eps)
{
  for (int i=2; i<=MAX_SUBTYPE; ++i)
    {
      if (fabs(pa.u[i] - prev_pa.u[i]) > eps)
        return false;
    }

  if (fabs(pa.beta[0][0] - prev_pa.beta[0][0]) > eps)
    return false;

  for (int i=1; i<=NONLEAF; ++i)
    for (int j=0; j<MAX_CHILD; ++j)
      if (fabs(pa.beta[i][j] - prev_pa.beta[i][j]) > eps)
        return false;

  for (int c=1; c<=MAX_COPY; ++c)
    if (fabs(pa.pi[c] - prev_pa.pi[c]) > eps)
      return false;

  for (int c=1; c<=MAX_COPY; ++c)
    for (int d=1; d<=c; ++d)
      if (fabs(pa.kappa[c][d] - prev_pa.kappa[c][d]) > eps)
        return false;

  for(int i=1; i<=MAX_SUBTYPE; ++i)
    if (fabs(pa.xi[0] - prev_pa.xi[i]) > eps)
      return false;

  for (int h=0; h<MAX_H; ++h)
    if (fabs(pa.omega[h] - prev_pa.omega[h]) > eps)
      return false;

  return true;
}

void write_params(std::ofstream& f, params pa)
{
  f << "u" << endl;
  for (int i=2; i<=MAX_SUBTYPE; ++i)
    f << pa.u[i] << "\t";
  f << endl;
  
  f << endl << "beta" << endl;
  pa.beta[0][0] = 0.5;
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
  f << endl;
  
  return;
}

int main(int argc, char* argv[]) {
  ifstream infile;
  infile.open("reads.txt");

  READS res;
  while(!infile.eof())
    {
      READ re;
      infile >> re.first >> re.second;
      res.push_back(re);
    }
  infile.close();
  
  Lbfgsb minimizer;
  minimizer.set_eps(1.0e-4);
  minimizer.set_maxit(30);
  vector<double> x0(2*MAX_SUBTYPE-1, 0.5);

  int em_maxit = 10;
  double em_eps = 1.0e-1;
  
  params pa_old; // pa_old initialization is needed
  init_params(pa_old);
  
  params pa_new;
  hyperparams hpa; // hpa initialization is needed
  init_hyperparams(hpa);
  
  for (int i=0; i<em_maxit; ++i)
    {
      maximization(pa_old, pa_new, hpa, res, minimizer); // res io is needed
      if (check_converge(pa_new, pa_old, em_eps)) // distance implementation is needed
        break;
      pa_old = pa_new;
    }

  ofstream outfile;
  outfile.open("optimized_params.txt");
  write_params(outfile, pa_new);
  outfile.close();
  
  return 0;
}
