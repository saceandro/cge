#include <fstream>
#include <math.h>
using namespace std;

double calc_dx_sigmoid(double x)
{
  if (x < -10.0)
    return exp(x) / (2.0*exp(x) + 1.0);
  else if (x < 10.0)
    return 1.0 / (exp(-x) + 2.0 + exp(x));
  else
    return exp(-x) / (2.0*exp(-x) + 1.0);
}

int main()
{
  ofstream ofile("sigmoid.txt");
  ofile << scientific;
  for (int i=-1000; i<=1000; ++i)
    {
      double x = ((double)i)/10;
      ofile << x << "\t" << calc_dx_sigmoid(x) << endl;
    }
  ofile.close();
}
