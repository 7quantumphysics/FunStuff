#ifndef RungeKutta_
#define RungeKutta_

#include <vector>
#include <functional> // Requires C++ 11

// To Compile:  g++ <filename>.cpp --std=c++11

/*
	Enter into RK: 
	    the function, initial function value,
	    derivative of function, initial value of derivative of function,
	    # of data points you wish to use, x_min, x_max, 
	    and the second derivative as a [C++] function of the independant variable 
	    and the function in question
*/

void RK(std::vector<double>& f, double& f_initial,
	std::vector<double>& df_dt, double& df_dt_initial,
	int& nDataPoints, double& min, double& max,
	const std::function <double (double, double)> SECOND_DERIVATIVE);

void RK2D(std::vector<double>& f1, double& f1_initial, std::vector<double>& f2, double& f2_initial,
	  std::vector<double>& df1_dt, double& df1_dt_initial, std::vector<double>& df2_dt, double& df2_dt_initial,
	  int& nDataPoints, double& min, double& max,
	  const std::function <double (double, double, double)> SECOND_DERIVATIVE1,
	  const std::function <double (double, double, double)> SECOND_DERIVATIVE2); 

#endif
