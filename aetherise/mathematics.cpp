
#include "mathematics.h"

#include <cmath>
#include <algorithm>

namespace aether {


Vector3 operator+(const Vector3& a, const Vector3& b)
{
	return {a.x+b.x, a.y+b.y, a.z+b.z};
}

Vector3 operator-(const Vector3& a, const Vector3& b)
{
	return {a.x-b.x, a.y-b.y, a.z-b.z};
}

Vector3 operator*(const Vector3& a, real s)
{
	return {a.x*s, a.y*s, a.z*s};
}

Vector3 operator/(const Vector3& a, real s)
{
	return {a.x/s, a.y/s, a.z/s};
}

real dot_product(const Vector3& a,const Vector3& b)
{
	return a.x*b.x + a.y*b.y + a.z*b.z;
}

real length(const Vector3& a)
{
	return std::sqrt(a.x*a.x + a.y*a.y + a.z*a.z);
}




Vector3 rotate_x(const Vector3& p,real w)
{
	using namespace std;
	
	return {	p.x,
				p.y*cos(w) - p.z*sin(w),
				p.z*cos(w) + p.y*sin(w)};
}

Vector3 rotate_y(const Vector3& p,real w)
{
	using namespace std;
	
	return {	p.x*cos(w) - p.z*sin(w),
				p.y,
				p.z*cos(w) + p.x*sin(w)};
}

Vector3 rotate_z(const Vector3& p,real w)
{
	using namespace std;
	
	return {	p.x*cos(w) - p.y*sin(w),
				p.y*cos(w) + p.x*sin(w),
				p.z};
}




real period_360(real deg)
{
	return norm_period(deg,real(360.0));
}



real period_2pi(real rad)
{
	return norm_period(rad,real(AETHER_2PI));
}



real period_24(real h)
{
	return norm_period(h,real(24.0));
}



bool approximate(real a, real b,real ratio)
{
	if (a==b)
		return true;

	if (a==0)
		return std::abs(b)<=ratio;
	if (b==0)
		return std::abs(a)<=ratio;
	auto d = std::abs(b-a);
	if (std::isinf(d))
		return false;
	return d <= ratio*std::min(std::abs(a),std::abs(b));

}


bool approximate_delta(real a, real b,real delta)
{
	return std::abs(a-b) <= std::abs(delta);
}

bool approximate_delta_periodic(real a, real b,real delta,real period)
{
	period = std::abs(period);
	delta = std::abs(delta);
	a = norm_period(a,period);
	b = norm_period(b,period);
	auto D = std::abs(a-b);
	return D <= delta || period-D <= delta;
}



real chi_square_pdf(real X,int n)
{
	auto n2 = n*0.5L;
	//return (std::pow(X,n2-1)*std::exp(-X*0.5)) / (std::pow(2,n2)*std::tgamma(n2));
	X *= 0.5L;
	return 0.5L*(std::pow(X,n2-1)/std::tgamma(n2)) * std::exp(-X);
}



real chi_square_cdf(real X,int n)
{
	// https://de.wikipedia.org/wiki/Chi-Quadrat-Verteilung

	const real x2 = X*0.5L;
	const int n2 = n/2;

	real c;
	if (n % 2) { // odd?
		real s = 0;
		for (int k=0;k<n2;k++) {
			s += std::pow(x2,k+0.5) / std::tgamma(k+1.5L);
		}
		c = std::erf(std::sqrt(x2)) - std::exp(-x2)*s;
	}
	else {
		real s = 0;
		for (int k=0;k<n2;k++) {
			s += std::pow(x2,k) / std::tgamma(k+1.0L);
		}
		c = 1 - std::exp(-x2)*s;
	}

	return c;
}



real integrate(real a,real b,std::function<real(real x)> f,int n,
			   Romberg_Cancelation_t cancel)
{
	// Lothar Papula: Mathematische Formelsammlung, 4. Auflage, S. 138

	if (n<0)
		n=0;
	if (n>31)
		n=31; // j is int

	std::vector<real> Tk(n);

	auto d = b-a;
	auto Tik = 0.5*d*(f(a)+f(b));

	real z = 0.5;
	for (int i=0;i<n;i++) {
		Tk[i] = Tik;
		if (cancel(i,Tk))
			break;

		KahanSum<real> s;
		z *= 2;
		auto d2z = d/(2*z);
		for (int j=1;j<=z;j++) {
			s += f(a + (2.0*j-1)*d2z);
		}

		Tik = 0.5*(Tk[0] + d/z*s.sum());

		real q = 1.0;
		for (int k=1;k<=i+1;k++) {
			q *= 4;
			auto tik = (q*Tik - Tk[k-1])/(q-1);
			Tk[k-1] = Tik;
			Tik = tik;
		}

		if (std::isinf(Tik) || std::isnan(Tik))
			break;
	}

	return Tik;
}




Polar3 polar3(const Vector3& v)
{
	Polar3 p;
	p.r = std::sqrt(sqr(v.x)+sqr(v.y)+sqr(v.z));
	p.th = std::atan2(std::hypot(v.x,v.y), v.z);
	p.ph = std::atan2(v.y, v.x);
	return p;
}



Vector3 vector3(const Polar3& p)
{
	Vector3 v;
	v.x = p.r*std::sin(p.th)*std::cos(p.ph);
	v.y = p.r*std::sin(p.th)*std::sin(p.ph);
	v.z = p.r*std::cos(p.th);
	return v;
}


real normal_pdf(real m,real s,real x)
{
	return 1/(std::sqrt(2*AETHER_PI_L)*s) * std::exp(-0.5L*sqr((x-m)/s));
}


real normal_cdf(real m,real s,real x)
{
	return 0.5*(1 + std::erf((1/std::sqrt(2)) * (x-m)/s));
}



double test_for_normality(const std::vector<double>& x,ADTestType transformation)
{
	// M. A. Stephens (1974): EDF Statistics for Goodness of Fit and Some Comparisons
	// Journal of the American Statistical Association, 69:347, 730-737

	size_t n = x.size();
	double m = mean_value(x.begin(),x.end());
	double s = sample_standard_deviation(x.begin(),x.end(),m);

	double A = 0;
	for (size_t i=1;i<=n;i++) {
		auto wi = (x[i-1]-m)/s;
		double zi = normal_cdf(0,1,wi);

		auto wk = (x[n-i]-m)/s;
		double zk = normal_cdf(0,1,wk);

		A += (2*i-1)*(std::log(zi) + std::log(1-zk));
	}

	A = -A/n - n;

	// Modification for A² if µ and σ are estimated.
	switch(transformation) {	
	case ADTestType::Stephens:
		// This transformation is most of the time > significance level.		
		// 10%=0.656, 5%=0.787, 1%=1.092
		A = A*(1. + (4. - 25./n)/n); // Stephens (1974)
		break;
	case ADTestType::DAgostino:
		// This transformation is more accurate, but half of the time below the significance level.
		// 10%=0.631, 5%=0.752, 1%=1.035
		A = A*(1. + (0.75 + 2.25/n)/n); // D'Agostino (1986)
		break;
	default:
		A = NAN;
	}
	return A;
}



double p_value_A(double A,ADTestType type)
{
	double p;

	if (type == ADTestType::DAgostino) {
		// Goodness-Of-Fit-Techniques, D'Agostino, Table 4.9, p. 127
		// ISBN: 0824774876

		if (A < 0.2) {
			p = 1-std::exp(-13.436 + (101.14 - 223.73*A)*A);
		}
		else if (A < 0.34) {
			p = 1-std::exp(-8.318 + (42.796 - 59.938*A)*A);
		}
		else if (A < 0.6) {
			p = std::exp(0.9177 + (-4.279 - 1.38*A)*A);
		}
		else {
			p = std::exp(1.2937 + (-5.709 - 0.0186*A)*A);
		}
	}
	else
		p = NAN;

	return p;
}



const ADTestQuantiles& test_quantiles(ADTestType type) {
	switch(type) {
	case ADTestType::Stephens:
		return ADTestStephensQuantiles;		
	case ADTestType::DAgostino:
		return ADTestDAgostinoQuantiles;		
	default:
		throw std::invalid_argument("unknown ADTest type");
	}
}



ConfidenceInterval confidence_interval(int k,int n)
{
	// https://de.wikipedia.org/wiki/Konfidenzintervall_f%C3%BCr_die_Erfolgswahrscheinlichkeit_der_Binomialverteilung
	const double c = 1.96; // 95% confidence

	auto k_= k+c*c*0.5;
	auto n_ = n + c*c;
	auto p_ = k_/n_;
	auto u_ = c*std::sqrt(p_*(1-p_)/n_);

	return {p_,u_};
}



double periodic_distance(double a,double b,double period)
{
	a = norm_period(a,period);
	b = norm_period(b,period);
	auto d = std::abs(b-a);
	return std::min(d,period-d);
}


}//aether
