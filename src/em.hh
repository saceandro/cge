#include <iostream>
#include <vector>
#include <utility>
#include <boost/math/distributions/beta.hpp>
#include <cmath>

#define FRACTIONS 10

typedef std::pair<int,int> read;
typedef std::vector<read> reads;

class state 
{
public:
  int k;
  int l[I];
  int r[I];
  int q;
  int h;
  double s;
  responsibility resp;
  
  state (int, int*, int*, int, int, double);
};

typedef std::vector<state> states;

class params
{
public:
  double u[I];
  double beta[I][J];
  double pi[C];
  double kappa[C][C];
  double xi[I];
  double omega[MAX_H];
  double beta_tilda[I];

  params (double*, double**, double*, double**, double*, double*, double*);
};
