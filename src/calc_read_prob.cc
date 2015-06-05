#include "calc_read_prob.hh"

double beta_pdf(double x, double a, double b)
{
  boost::math::beta_distribution<> beta(a, b);
  return pdf(beta, x);
}

double prob_u(std::vector<double> u, std::vector<double> a, std::vector<double>b)
{
  double p = 1;
  for (std::vector<double>::iterator i = u.begin(), j = a.begin(), k = b.begin(); i != u.end(); ++i, ++j, ++k)
    {
      p *= beta_pdf(*i, *j, *k);
    }
  return p;
}

int main()
{
  double u_i, a_u, b_u;
  std::cin >> u_i >> a_u >> b_u;

  std::cout << beta_pdf(u_i, a_u, b_u) << std::endl;
  return 0;
}
