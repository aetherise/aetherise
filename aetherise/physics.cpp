#include "physics.h"

#include "mathematics.h"

#include <cmath>

namespace aether {


namespace {

double compressibility_of_air(double p, double T,double xw)
{
	// Philip E. Ciddor: "Refractive index of air: new equations for the visible and near infrared", 1996
	// referencing
	// R. S. Davis, "Equation for the determination of the density of moist air", (1981/1991),
	// Metrologia 29, 67–70 (1992)

	const auto a0 {1.58123e-6}; // (K/Pa)
	const auto a1 {-2.9331e-8}; // (1/Pa)
	const auto a2 {1.1043e-10}; // (1/(KPa))
	const auto b0 {5.707e-6}; // (K/Pa)
	const auto b1 {-2.051e-8}; // (1/Pa)
	const auto c0 {1.9898e-4}; // (K/Pa)
	const auto c1 {-2.376e-6}; // (1/Pa)
	const auto d {1.83e-11}; // (K²/Pa²)
	const auto e {-0.765e-8}; // (K²/Pa²)

	auto t = in_Celsius(T);
	auto pT = p/T;
	return 1-pT*((a0 + (a1 + a2*t)*t + ((b0+b1*t) + (c0+c1*t)*xw)*xw) - pT*(d+e*xw*xw));
}

}


double refractive_index_of_air(double p,double T,double lambda,double h,double xc)
{
	// Philip E. Ciddor: "Refractive index of air: new equations for the visible and near infrared", 1996

	const auto A {1.2378847e-5}; // (1/K²)
	const auto B {-1.9121316e-2}; // (1/K)
	const auto C {33.93711047};
	const auto D {-6.3431645e3}; // (K)

	// svp: saturation vapor pressure of water vapor
	auto svp = std::exp((A*T + B)*T + C + D/T); // (Pa)

	const auto a {1.00062};
	const auto b {3.14e-8}; // (1/Pa)
	const auto g {5.6e-7}; // (1/°C²)


	/*auto fxw = [&](double p,double t,double h){
		auto f = a + b*p + g*t*t;
		auto xw = f*h*svp/p;
		return xw;
	};*/


	auto t = in_Celsius(T); // K to °C	
	auto f = a + b*p + g*t*t; // enhancement factor of water vapor in air
	auto xw = f*h*svp/p; // molar fraction of water vapor in moist air


	const auto k0 {238.0185}; // (1/µm²)
	const auto k1 {5792105.}; // (1/µm²)
	const auto k2 {57.362}; // (1/µm²)
	const auto k3 {167917.}; // (1/µm²)
	auto s2 = sqr(1./lambda); // (1/µm)
	auto n_as1 = (k1/(k0-s2) + k3/(k2-s2))*1e-8;
	auto n_axs = 1 + n_as1*(1 + 0.534e-6*(xc-450));

	const auto w0 {295.235}; // (1/µm²)
	const auto w1 {2.6422}; // (1/µm²)
	const auto w2 {-0.032380}; // (1/µm^4)
	const auto w3 {0.004028}; // (1/µm^6)
	const auto cf {1.022};
	auto n_ws = 1 + cf*(w0 + s2*(w1 + s2*(w2 + w3*s2)))*1e-8;

	// molar mass of dry air containing xc ppm of CO2
	auto M_a = 0.0289635 + 12.011e-9*(xc-400); // (kg/mol)

	// compressibility of dry air
	//auto Z_a = compressibility_of_air(101325,288.15,0);
	auto Z_a = 0.999592211536;

	// compressibility of pure water vapor
	//auto Z_w = compressibility_of_air(1333,293.15,1);

	// M_w: molar mass of water vapor
	const auto M_w {0.018015}; // (kg/mol)
	const double R {Gas_Constant}; // (J/(mol K))

	auto r_axs = 101325./288.15*M_a/(R*Z_a);
	//auto r_ws = 1333./293.15*M_w/(R*Z_w); // wrong?
	//auto r_ws = 0.00985938; // value from NIST web page
	auto r_ws = 0.009859437276649; // calculated with Z_w

	auto Z = compressibility_of_air(p,T,xw);

	auto r_a = p*M_a*(1-xw)/(Z*R*T);
	auto r_w = p*M_w*   xw /(Z*R*T);

	auto L_a = (sqr(n_axs)-1) / (sqr(n_axs)+2);
	auto L_w = (sqr(n_ws) -1) / (sqr(n_ws) +2);
	auto L = (r_a/r_axs)*L_a + (r_w/r_ws)*L_w;
	return std::sqrt((1+2*L) / (1-L)); // second method
	//return 1+(r_a/r_axs)*(n_axs-1) + (r_w/r_ws)*(n_ws-1); // first method
}



double barometric_formula(double h,double p0,double T0)
{
	// https://de.wikipedia.org/wiki/Barometrische_H%C3%B6henformel
	return p0 * std::pow(1-0.0065*h/T0, 5.255);
}


}//
