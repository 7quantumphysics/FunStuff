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

void RK2D(vector<double>& f1, double& f1_initial, vector<double>& f2, double& f2_initial,
	  vector<double>& df1_dt, double& df1_dt_initial, vector<double>& df2_dt, double& df2_dt_initial,
	  int& nDataPoints, double& min, double& max,
	  const function <double (double, double, double)> SECOND_DERIVATIVE1,
	  const function <double (double, double, double)> SECOND_DERIVATIVE2)
{
  double t = min; // initialize independant variable 
  double dt = (max-min)/nDataPoints; // Define the small interval
  f1.push_back(f1_initial); // Initialize 1st dependant variable
  f2.push_back(f2_initial); // Initialize 2nd dependant variable
  df1_dt.push_back(df1_dt_initial); // Initialize first derivative of 1st dependant variable
  df2_dt.push_back(df2_dt_initial); // Initialize first derivative of 2nd dependant variable
  
  for(unsigned int i = 0; t <= max; i++)
    {
      double k1f1 = (df1_dt[i])*dt;
      double k1f2 = (df2_dt[i])*dt;
      double k1g1 = SECOND_DERIVATIVE1(t, f1[i], f2[i])*dt;
      double k1g2 = SECOND_DERIVATIVE2(t, f1[i], f2[i])*dt;

      double k2f1 = (df1_dt[i] + 0.5*k1g1)*dt;
      double k2f2 = (df2_dt[i] + 0.5*k1g2)*dt;
      double k2g1 = SECOND_DERIVATIVE1(t+0.5*dt, f1[i]+0.5*k1f1, f2[i]+0.5*k1f2)*dt;
      double k2g2 = SECOND_DERIVATIVE2(t+0.5*dt, f1[i]+0.5*k1f1, f2[i]+0.5*k1f2)*dt;

      double k3f1 = (df1_dt[i] + 0.5*k2g1)*dt;
      double k3f2 = (df2_dt[i] + 0.5*k2g2)*dt;
      double k3g1 = SECOND_DERIVATIVE1(t+0.5*dt, f1[i]+0.5*k2f1, f2[i]+0.5*k2f2)*dt;
      double k3g2 = SECOND_DERIVATIVE2(t+0.5*dt, f1[i]+0.5*k2f1, f2[i]+0.5*k2f2)*dt;

      double k4f1 = (df1_dt[i] + 0.5*k3g1)*dt;
      double k4f2 = (df2_dt[i] + 0.5*k3g2)*dt;
      double k4g1 = SECOND_DERIVATIVE1(t+dt, f1[i]+k3f1, f2[i]+k3f2)*dt;
      double k4g2 = SECOND_DERIVATIVE2(t+dt, f1[i]+k3f1, f2[i]+k3f2)*dt;

      // Declare incremental step for f and its derivative using k's //

      double step_f1 = (1.0/6.0)*(k1f1 + 2.0*k2f1 + 2.0*k3f1 + k4f1);
      double step_f2 = (1.0/6.0)*(k1f2 + 2.0*k2f2 + 2.0*k3f2 + k4f2);
      double step_g1 = (1.0/6.0)*(k1g1 + 2.0*k2g1 + 2.0*k3g1 + k4g1);
      double step_g2 = (1.0/6.0)*(k1g2 + 2.0*k2g2 + 2.0*k3g2 + k4g2);

      // Update the function, its derivative and the independant variable //

      f1.push_back(f1[i]+step_f1);
      f2.push_back(f2[i]+step_f2);
      df1_dt.push_back(df1_dt[i]+step_g1);
      df2_dt.push_back(df2_dt[i]+step_g2);
      t += dt;
    }
}
