#include "RungeKutta.h"
#include <iostream>
#include <cmath>
#include <functional>
#include <vector>

using namespace std;

/* -------------------------------------------------------- */
/* //    Define Functions' Second Order Derivatives      // */
/* //    in terms of the Independant Variable and the    // */
/* //    Function itself                                 // */
/* -------------------------------------------------------- */
/* i.e.)                                                    */

double d2x_dt2(double t, double x, double y)
{
  double denom = pow( (x*x + y*y),(3./2.) );
  double x_double_dot = (-4.*M_PI*M_PI/denom)*(x);
  return x_double_dot;
}

double d2y_dt2(double t, double x, double y)
{
  double denom = pow( (x*x + y*y),(3./2.) );
  double y_double_dot = (-4.*M_PI*M_PI/denom)*(y);
  return y_double_dot;
}

int main()
{
  /* ---------------------------------------------------------- */
  /* //    Declare the Functions, First Derivatives, and     // */
  /* //           their Respective Initial Values            // */
  /* ---------------------------------------------------------- */
  /* i.e.)                                                      */

  vector<double> X;
  vector<double> X_dot;
  vector<double> Y;
  vector<double> Y_dot;
  
  double X_initial = 1.0;
  double X_dot_initial = 0.0;
  double Y_initial = 0.0;
  double Y_dot_initial = 2.*M_PI;

  /* ------------------------------------------------------------- */
  /* //    Declare the Number of Data Points and the Minimum    // */
  /* //     and Maximum Values of the Independant Variable      // */
  /* ------------------------------------------------------------- */
  /* i.e.)                                                         */

  int nPoints = 100000;
  double minimum = 0.0;
  double maximum = 10.0;

  /* -------------------------------------------------------------- */
  /* //     Re-Declare the Second Derivative as a Functional     // */
  /* -------------------------------------------------------------- */
  /* i.e.)                                                          */

  function <double (double, double, double)> Second_Derivative1 = d2x_dt2;
  function <double (double, double, double)> Second_Derivative2 = d2y_dt2;

  /* ------------------------------------------ */
  /* //     Call the RungeKutta Function     // */
  /* ------------------------------------------ */
  /* i.e.)                                      */

  RK(X, X_initial, X_dot, X_dot_initial,
     Y, Y_initial, Y_dot, Y_dot_initial,
     nPoints, minimum, maximum,
     Second_Derivative1, Second_Derivative2);

  /* -------------------------------------------------------- */
  /* //      Perform analysis here using the desired       // */
  /* //      information from the RungeKutta function      // */
  /* -------------------------------------------------------- */
    
  return 0;
}
