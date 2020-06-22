#include "RungeKutta.h"
#include <iostream>
#include <cmath>
#include <functional>
#include <vector>

using namespace std;

/* ----------------------------------------------_--------- */
/* //    Define a Function's Second Order Derivative     // */
/* //    in terms of the Independant Variable and the    // */
/* //    Function itself                                 // */
/* -------------------------------------------------------- */
/* i.e.)                                                    */

double dt2_dt2(double t, double y)
{
  double c = 0.2;
  double F = 0.5;
  double W = 0.7;
  double y_double_dot = -c*sin(y) + F*cos(W*t);
  return y_double_dot;
}

int main()
{
  /* -------------------------------------------------------- */
  /* //    Declare the Function, First Derivative, and     // */
  /* //          their Respective Initial Values           // */
  /* -------------------------------------------------------- */
  /* i.e.)                                                    */
  
  vector<double> Y;
  vector<double> Y_dot;
  double Y_initial = 0.1;
  double Y_dot_initial = 0.0;

  /* ------------------------------------------------------------- */
  /* //    Declare the Number of Data Points and the Minimum    // */
  /* //     and Maximum Values of the Independant Variable      // */
  /* ------------------------------------------------------------- */
  /* i.e.)                                                         */

  int nPoints = 100000;
  double minimum = 0.0;
  double maximum = 200.0;

  /* -------------------------------------------------------------- */
  /* //     Re-Declare the Second Derivative as a Functional     // */
  /* -------------------------------------------------------------- */
  /* i.e.)                                                          */

  function <double (double, double)> Second_Derivative = dt2_dt2;

  /* ------------------------------------------ */
  /* //     Call the RungeKutta Function     // */
  /* ------------------------------------------ */
  /* i.e.)                                      */

  RK(Y, Y_initial, Y_dot, Y_dot_initial,
     nPoints, minimum, maximum, Second_Derivative);

  /* -------------------------------------------------------- */
  /* //      Perform analysis here using the desired       // */
  /* //      information from the RungeKutta function      // */
  /* -------------------------------------------------------- */
    
  return 0;
}
