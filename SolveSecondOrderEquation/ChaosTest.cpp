#include "RungeKutta.h"
#include <iostream>
#include <cmath>
#include <functional>
#include <vector>

using namespace std;

double d2f_dt2(double time, double function)
{
  double c = 0.2;
  double F = 0.5;
  double W = 0.7;
  double ans = -c*sin( function ) + F*cos( W*time );

  return ans;
}

int main()
{
  vector<double> theta; // Theta as a function of time; the thing that evolves with time
  double theta_initial = 0.1; // Initial value of theta
  vector<double> omega; // Omega; First derivative of theta with respect to time written as its own vector variable
  double omega_initial = 0.0; // Initial value of omega
  double duration = 100000.0; // The higher this is, the more precise the data
  double minimum = 0.0; // Beginning of x-axis
  double maximum = 200.0; // End of x-axis

  function <double (double, double)> SD = d2f_dt2; // Declared only to match the definition in the header file
  
  RK(theta, theta_initial,
     omega, omega_initial,
     duration, minimum, maximum,
     SD);
  
  return 0;
}
