#include "setting.hh"
#include "em.hh"

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

double calc_alpha(double t_i, double n_i, int N)
{
  return log(N * n_i) / t_i;
}

double read_prob(read re, state st, double* n, double epsilon)
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
          int parent = (i + J - 2) / J;
          while (parent > st.q)
            {
              prev_parent = parent;
              parent = (prev_parent + J - 2) / J;
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
      
      numerator += epsilon * n[i] * ((1 - x) * st.l[i] + x * (st.l[i] - st.r[i])) + (1 - epsilon) * n[i] * x * st.r[i];
    }
  
  double mu = numerator / denominator;

  boost::math::binomial_distribution<> binom(re.first, mu);

  return pdf(binom, re.second);
}

double variant_fraction(double x_i, int h, int i, double n_i, double t_i, double beta_tilda)
{
  return n_i / x_i / beta_tilda;
}

double variant_fraction_alpha(double x_i, int h, int i, double alpha_i, double t_i, double beta_tilda, double N)
{
  return exp(alpha_i * t_i)  / x_i / beta_tilda / N;
}

double d_alpha_variant_fraction(double x_i, int h, int i, double alpha_i, double t_i, double beta)
{
  return t_i * variant_fraction(x_i, h, i, alpha_i, t_i, beta);
}

double d_t_variant_fraction(double x_i, int h, int i, double alpha_i, double t_i, double beta)
{
  return alpha_i * variant_fraction(x_i, h, i, alpha_i, t_i, beta);
}

double responsibility_numerator(state st, params pa_old, read *re)
{
  double product = 1;
  for (int i=1; i<I; i++)
    product *= pa_old.pi[st.l[i]] * pa_old.kappa[st.l[i]][st.r[i]];
  product *= pa_old.xi[st.q] * pa_old.omega[st.h] * variant_fraction(st.s, st.h, st.q, calc_n(st.q, pa_old.beta), calc_t(st.q, pa_old.u), pa_old.beta_tilda[st.q]) * read_prob(re, st, calc_n(st.q, pa_old.beta), epsilon);
  return product;
}


double responsibility_partition(int k, int K, params pa_old, read *re)
{
  state curr_st;
  
  double e = 0;

  curr_st.k = k;
  curr_st.l[0] = 2;
  curr_st.r[0] = 0;
  double d[5];
  double c[5];
      
  for (curr_st.l[1]=0; curr_st.l[1]<C; ++curr_st.l[1])
    {
      double d[1] = 0;

      for (curr_st.r[1]=0; curr_st.r[1]<curr_st.l[1]; ++curr_st.r[1])
        {
          double c[1] = 0;

          for (curr_st.l[2]=0; curr_st.l[2]<C; ++curr_st.l[2])
            {
              double d[2] = 0;

              for (curr_st.r[2]=0; curr_st.r[2]<curr_st.l[2]; ++curr_st.r[2])
                {
                  double c[2] = 0;

                  for (curr_st.l[3]=0; curr_st.l[3]<C; ++curr_st.l[3])
                    {
                      double d[3] = 0;

                      for (curr_st.r[3]=0; curr_st.r[3]<curr_st.l[3]; ++curr_st.r[3])
                        {
                          double c[3] = 0;
                      
                          for (curr_st.l[4]=0; curr_st.l[4]<C; ++curr_st.l[4])
                            {
                              double d[4] = 0;

                              for (curr_st.r[4]=0; curr_st.r[4]<curr_st.l[4]; ++curr_st.r[4])
                                {
                                  double c[4] = 0;
                      
                                  for (curr_st.q=0; curr_st.q<I; ++curr_st.q)
                                    {
                                      double b = 0;
                      
                                      for (curr_st.h=0; curr_st.h<2**J; ++curr_st.h)
                                        {
                                          double a = 0;
                          
                                          for (int s=0; s<FRACTIONS; ++s)
                                            {
                                              curr_st.s = ((double)s) / FRACTIONS;
                                              a += variant_fraction(curr_st.s, curr_st.h, curr_st.q, calc_n(curr_st.q, pa_old.beta), calc_t(curr_st.q, pa_old.u), pa_old.beta_tilda[curr_st.q]) * read_prob(re, curr_st, calc_n(curr_st.q, pa_old.beta), epsilon);
                                            }
                                          b += pa_old.omega[h] * a;
                                        }
                                      c[4] += pa_old.xi[q] * b;
                                    }
                                  d[4] += pa_old.kappa[curr_st.l[4]][curr_st.r[4]] * c[4];
                                }
                              c[3] += pa_old.pi[curr_st.l[4]] * d[4];
                            }
                          d[3] += pa_old.kappa[curr_st.l[3]][curr_st.r[3]] * c[3];
                        }
                      c[2] += pa_old.pi[curr_st.l[3]] * d[3];
                    }
                  d[2] += pa_old.kappa[curr_st.l[2]][curr_st.r[2]] * c[2];
                }
              c[1] += pa_old.pi[curr_st.l[2]] * d[2];
            }
          d[1] += pa_old.kappa[curr_st.l[1]][curr_st.r[1]] * c[1];
        }
      e += pa_old.pi[curr_st.l[1]] * d[1];
    }
  return e;
}

void responsibility(states sts, params old_pa, int K, reads res)
{
  state st = new state;
  for (int k=0, std::vector<read>::iterator re = res.begin(); k<K; ++k, ++re)
    {
      double resp_part = responsibility_partition(k, K, old_pa, re);
      int l[5];
      int r[5];
      
      for (l[1]=0; l[1]<C; ++l[1])
        {
          for (r[1]=0; r[1]<l[1]; ++r[1])
            {
              for (l[2]=0; l[2]<C; ++l[2])
                {
                  for (r[2]=0; r[2]<l[2]; ++r[2])
                    {
                      for (l[3]=0; l[3]<C; ++l[3])
                        {
                          for (r[3]=0; r[3]<l[3]; ++r[3])
                            {
                              for (l[4]=0; l[4]<C; ++l[4])
                                {
                                  for (r[4]=0; r[4]<l[4]; ++r[4])
                                    {
                                      for (int q=0; q<I; ++q)
                                        {
                                          for (int h=0; h<2**J; ++h)
                                            {
                                              for (int s=0; s<FRACTIONS; ++s)
                                                {
                                                  state st = new state;
                                                  
                                                  st.k = k;
                                                  
                                                  for (int i=1; i<=I; ++i)
                                                    {
                                                      st.l[i] = l[i];
                                                      st.r[i] = r[i];
                                                    }

                                                  st.q = q;
                                                  st.h = h;
                                                  st.s = ((double)s) / FRACTIONS;

                                                  st.resp = responsibility_numerator(st, old_pa, re) / resp_part;
                                                  sts.push_back(st);
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void maximization(params old_pa, params new_pa, int K, reads res)
{
  states sts = new states;
  responsibility(sts, old_pa, K, res);

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

}
