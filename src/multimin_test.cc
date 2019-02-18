#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
using namespace std;

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h> // gsl vector stuff
#include <gsl/gsl_multimin.h> // gsl multidimensional minimization
// function prototypes
inline double sqr (double x) 
{
  return x*x;
}
;
// square a double
double my_f (const gsl_vector *xvec_ptr, void *params);

void my_df (const gsl_vector *xvec_ptr, void *params, gsl_vector *df_ptr);

void my_fdf (const gsl_vector *xvec_ptr, void *params_ptr, double *f_ptr, gsl_vector *df_ptr);

//*********************************************************************
int
main ()
{

  // Exact position of the twoédimensional minimum (1,2).
  double minima[2] = 
    {
      1.0, 2.0 
    }
  ;

  // define and set up the gsl multimin function
  gsl_multimin_function_fdf my_func;

  my_func.f = &my_f;
  // the function to be minimized
  my_func.df = &my_df;
  // the gradient of the function
  my_func.fdf = &my_fdf;
  // combined (see the function below)
  my_func.n = 2;
  // size of x vectors
  my_func.params = &minima;
  // parameters available to the function
  // Allocate x vector and set starting point, e.g., x = (5,7)
  gsl_vector *xvec_ptr = gsl_vector_alloc (2);

  gsl_vector_set (xvec_ptr, 0, 5.0);

  gsl_vector_set (xvec_ptr, 1, 7.0);

  // allocate and set the minimizer and its type (see gsl manual)
   const gsl_multimin_fdfminimizer_type *type_ptr
     = gsl_multimin_fdfminimizer_conjugate_fr;

    gsl_multimin_fdfminimizer *minimizer_ptr
      = gsl_multimin_fdfminimizer_alloc (type_ptr, 2);

    // set the tolerance and starting step size
    double step_size = 0.01;
    double tolerance = 1.e-4;

    gsl_multimin_fdfminimizer_set (minimizer_ptr, &my_func, xvec_ptr,
                                   step_size, tolerance);

    size_t iteration = 0;
    // initialize iteration counter
    size_t max_iterations = 100;
    // stop at this iteration if not converged
    int status = 0;
    cout << "iter x y value " << endl;
     do
       {

         iteration++;

         status = gsl_multimin_fdfminimizer_iterate (minimizer_ptr);

         if (status)
           {

             break;
             // this should only happen if an error code is returned
           }

         // check for convergence (state is either GSL_CONTINUE or GSL_SUCCESS)
         status = gsl_multimin_test_gradient (minimizer_ptr->gradient, tolerance);

         if (status == GSL_SUCCESS) // if we’re done, print out the details
           {

             cout << "Minimum found at: " << endl;

           }

         cout << setw(3) << iteration << " "
              << fixed << setprecision(5)
              << setw(10) << gsl_vector_get (minimizer_ptr->x, 0) << " "
              << setw(10) << gsl_vector_get (minimizer_ptr->x, 1) << " "
              << setw(12) << minimizer_ptr->f
              << endl;

       }

     while (status == GSL_CONTINUE && iteration < max_iterations);

     gsl_multimin_fdfminimizer_free (minimizer_ptr);
     // free the minimizer
     gsl_vector_free (xvec_ptr);
     // free the vector
     return 0;

}

//*********************************************************************
// Simple paraboloid centered on (dp[0],dp[1])
double
my_f (const gsl_vector *xvec_ptr, void *params)
{

  cout << "xvec_ptr: " << xvec_ptr->data[0] << "\t" << xvec_ptr->data[1] << endl;
  
  double *dp = (double *)params;


  double x = gsl_vector_get(xvec_ptr, 0);
  double y = gsl_vector_get(xvec_ptr, 1);

  return ( 10.0 * sqr(x - dp[0]) + 20.0 * sqr(y - dp[1]) + 30.0 );

}

// The gradient of f, df = (df/dx, df/dy).
void
my_df (const gsl_vector *xvec_ptr, void *params,
       gsl_vector *df_ptr)
{

  double *dp = (double *)params;

  double x = gsl_vector_get(xvec_ptr, 0);

  double y = gsl_vector_get(xvec_ptr, 1);

  gsl_vector_set(df_ptr, 0, 20.0 * (x - dp[0]));

  gsl_vector_set(df_ptr, 1, 40.0 * (y - dp[1]));

}

// Compute both f and df together.
void
my_fdf (const gsl_vector *xvec_ptr, void *params_ptr,
        double *f_ptr, gsl_vector *df_ptr)
{

  *f_ptr = my_f(xvec_ptr, params_ptr);

  my_df(xvec_ptr, params_ptr, df_ptr);

}

//*****************************************
