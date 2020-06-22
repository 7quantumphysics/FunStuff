#include "RungeKutta.h"
#include <functional> // Requires C++ 11
#include <iostream>
#include <cmath>
#include <vector>

// To Compile:  g++ <filename>.cpp --std=c++11

using namespace std;

void RK(vector<double>& f, double& f_initial,
	vector<double>& df_dt, double& df_dt_initial,
	int& nDataPoints, double& min, double& max,
       	const function <double (double, double)> SECOND_DERIVATIVE)
{
  double t = min; // initialize independant variable 
  double dt = (max-min)/nDataPoints; // Define the small interval
  f.push_back(f_initial); // Initialize dependant variable
  df_dt.push_back(df_dt_initial); // Initialize first derivative of dependant variable
  
  for(unsigned int i = 0; t <= max; i++)
    {
      double k1f = (df_dt[i])*dt;
      double k1g = SECOND_DERIVATIVE(t, f[i])*dt;

      double k2f = (df_dt[i] + 0.5*k1g)*dt;
      double k2g = SECOND_DERIVATIVE(t+0.5*dt, f[i]+0.5*k1f)*dt;

      double k3f = (df_dt[i] + 0.5*k2g)*dt;
      double k3g = SECOND_DERIVATIVE(t+0.5*dt, f[i]+0.5*k2f)*dt;

      double k4f = (df_dt[i] + 0.5*k3g)*dt;
      double k4g = SECOND_DERIVATIVE(t+dt, f[i]+k3f)*dt;

      // Declare incremental step for f and its derivative using k's //

      double step_f = (1.0/6.0)*(k1f + 2.0*k2f + 2.0*k3f + k4f);
      double step_g = (1.0/6.0)*(k1g + 2.0*k2g + 2.0*k3g + k4g);

      // Update the function, its derivative and the independant variable //

      f.push_back(f[i]+step_f);
      df_dt.push_back(df_dt[i]+step_g);
      t += dt;
    }
}

