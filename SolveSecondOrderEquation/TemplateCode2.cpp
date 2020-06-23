#include "RungeKutta.h"
#include <iostream>
#include <cmath>
#include <functional>
#include <vector>
#include <fstream>

using namespace std;
ofstream OrbitFile("Orbit.dat"); // Path of the earth
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
  vector<double> VX;
  vector<double> Y;
  vector<double> VY;
  
  double X_initial = 1.0; // Distance in A.U.
  double VX_initial = 0.0;
  double Y_initial = 0.0;
  double VY_initial = 2.*M_PI;

  /* ------------------------------------------------------------- */
  /* //    Declare the Number of Data Points and the Minimum    // */
  /* //     and Maximum Values of the Independant Variable      // */
  /* ------------------------------------------------------------- */
  /* i.e.)                                                         */

  int nPoints = 100000;
  double t_min = 0.0;
  double t_max = 1.0; // Time in Earth Years

  /* -------------------------------------------------------------- */
  /* //     Re-Declare the Second Derivative as a Functional     // */
  /* -------------------------------------------------------------- */
  /* i.e.)                                                          */

  function <double (double, double, double)> AX = d2x_dt2;
  function <double (double, double, double)> AY = d2y_dt2;

  /* ------------------------------------------ */
  /* //     Call the RungeKutta Function     // */
  /* ------------------------------------------ */
  /* i.e.)                                      */

  RK2D(X, X_initial, Y, Y_initial,
       VX, VX_initial, VY, VY_initial,
       nPoints, t_min, t_max,
       AX, AY);

  /* -------------------------------------------------------- */
  /* //      Perform analysis here using the desired       // */
  /* //      information from the RungeKutta function      // */
  /* -------------------------------------------------------- */
  
  for(int i = 0; i < nPoints; i++)
    {
      //cout << X[i] << "    " << Y[i] << endl;
      OrbitFile << X[i] << "    " << Y[i] << endl;
    }
    
  return 0;
}
