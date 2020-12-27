
#include "test.h"
#include "../mathematics.h"

#include <cmath>


namespace test {

std::vector<_Test> _tests;


bool approximate(double a, double b,double f)
{
	return aether::approximate(a,b,f);

}

bool approximate_delta(double a, double b,double delta)
{
	return aether::approximate_delta(a,b,delta);

}

bool approximate_delta_periodic(double a, double b,double delta,double period)
{
	return aether::approximate_delta_periodic(a,b,delta,period);

}

}//aether

