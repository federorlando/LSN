#include "TSP_position.h"

#include <cmath>
#define _USE_MATH_DEFINES



TSP_position::TSP_position()
{
	m_x = 0;
	m_x = 0;
}



TSP_position::TSP_position(double x, double y)
{
	m_x=x;
	m_y=y;
}



TSP_position::~TSP_position(){}



/*
double TSP_position::Eval(double x) const
{
	double alpha = 1./(m_sigma*sqrt(2.*M_PI));
	double power = -pow(x-m_mean,2)/(2.*pow(m_sigma,2));	
	return alpha*exp(power);
}
*/
