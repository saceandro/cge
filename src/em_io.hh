// #define BOOST_MATH_DOMAIN_ERROR_POLICY ignore_error
#include "setting.hh"
#include <iostream>
#include <fstream>
#include <vector>
#include <utility>
// #include <boost/math/distributions/beta.hpp>
// #include <boost/math/distributions/binomial.hpp>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gegenbauer.h>
#include <gsl/gsl_multimin.h>

#define FRACTIONS 10
#define GEGEN_MAX 200

#define calc_parent(i) (((i) + MAX_CHILD - 2)/MAX_CHILD)
#define calc_remainder(i) (((i) + MAX_CHILD - 2) % MAX_CHILD)
#define calc_child(i,reminder) (MAX_CHILD * ((i)-1) + (reminder) + 2)

typedef std::pair<int,int> READ;
typedef std::vector< READ* > READS;
typedef std::vector< double > V;
typedef std::vector< V > VV;

class state 
{
public:
  int k;
  std::vector<int> l;
  std::vector<int> r;
  int q;
  int h;
  double s;
  double resp;

  // use default constructor
  // state (int, std::vector<int>, std::vector<int>, int, int, double);
};

typedef std::vector<state> states;

class params
{
public:
  std::vector<double> u;
  std::vector<std::vector<double> > beta;
  std::vector<double> pi;
  std::vector<std::vector<double> > kappa;
  std::vector<double> xi;
  std::vector<double> omega;
  std::vector<double> beta_tilda;

  // use default constructor
  // params (std::vector<double>, std::vector<std::vector<double> >, double*, double**, std::vector<double>, std::vector<double>, std::vector<double>);
};

class hyperparams
{
public:
  double au;
  double bu;
  double abeta;
  double bbeta;
  double epsilon;
}
  ;
