/**
 * \file
 *
 * \~german
 * @brief Mathematische Typen und Funktionen
 *
 * \~english
 * @brief Mathematic types and functions
 *
 */

#ifndef MATHEMATICS_H
#define MATHEMATICS_H

#include <cmath>
#include <numeric>
#include <stdexcept>
#include <vector>
#include <functional>
#include <limits>

namespace aether {

using real = long double;





// 2π
#define AETHER_2PI		6.283185307179586477
// π
#define AETHER_PI		3.14159265358979323846
// π/2
#define AETHER_PI_2		1.57079632679489661923



// 2π
#define AETHER_2PI_L	6.283185307179586477L
// π
#define AETHER_PI_L		3.14159265358979323846L
// π/2
#define AETHER_PI_2_L	1.57079632679489661923L


constexpr real rad(real deg) {
	return AETHER_PI_L/180.0 * deg;
}

constexpr real deg(real rad) {
	return 180.0/AETHER_PI_L * rad;
}

constexpr real rad_to_h(real rad) {
	return 12.0/AETHER_PI_L*rad;
}

constexpr real h_to_rad(real h) {
	return AETHER_PI_L/12.0*h;
}



/**
 * \~german Mittelwert
 * Arithmetisches Mittel
 *
 * \~english
 * Arithmetic mean
 * 
 * \~
 * 
 */
template<typename Iter,typename T = typename std::iterator_traits<Iter>::value_type>
T mean_value(Iter it,Iter end)
{	
	static_assert(std::is_floating_point<T>::value,"floating point needed");

	auto n = std::distance(it,end);
	if (n <= 0)
		return T(NAN);

	T mean = *it;
	for (++it; it != end; ++it)
		mean += *it;
	return mean/n; // T/unsigned !
}



template<typename T,size_t N>
T mean_value(const std::array<T,N>& a)
{
	return mean_value(a.begin(),a.end());
}


/**
 * \~german Median
 * 
 * Bei gerader Anzahl an Werten wird der Mittelwert der beiden Werte 
 * in der Mitte berechnet.
 * 
 * \~english
 * If the number of values is even, then the mean value of the two
 * values in the middle is calculated.
 * 
 * \~
 * @param sorted_values must be sorted
 * @return
 */
template<typename T>
T median_value(const std::vector<T>& sorted_values)
{
	static_assert(std::is_floating_point<T>::value,"floating point needed");
	
	if (sorted_values.empty())
		return T(NAN);
	
	T m;	
	auto i = sorted_values.size()/2;
	if (sorted_values.size()%2) {
		m = sorted_values.at(i);	
	}
	else {				
		m = 0.5*(sorted_values.at(i-1)+sorted_values.at(i));		
	}
	return m;
}




/**
 * \~german
 * Kahan-Babuška-Verfahren
 *
 * Rundungsfehler werden berechnet und getrennt summiert.
 *
 * \~english
 * Kahan summation algorithm
 *
 * (compensated summation)
 */
template<typename T>
class KahanSum
{
	static_assert(std::is_floating_point<T>::value,"floating point needed");
private:
	T s;
	T w;

public:
	KahanSum()
		:s(),w()
	{
	}

	KahanSum(const T& s0)
		: s{s0},w()
	{
	}

	void add(const T& a) {
		T sm = s + a;
		w = w+(a+(s-sm));
		s = sm;
	}

	T sum() const {
		return s + w;
	}

	KahanSum<T>& operator+=(const T& a)
	{
		add(a);
		return *this;
	}

	// cast returns sum
	explicit operator T() const {
		return sum();
	}

};





//-------------
// Vektor
//-------------



/**
 * \~german
 * Euklidischer 3D-Vektor
 *
 * \~english
 * Euclidian 3D vector
 */
struct Vector3
{
	real x,y,z;
};

Vector3 operator+(const Vector3& a, const Vector3& b);

Vector3 operator-(const Vector3& a, const Vector3& b);

Vector3 operator*(const Vector3& a, real s);

Vector3 operator/(const Vector3& a, real s);



/**
 * \~german Skalarprodukt
 * \~english
 * \~
 * @param a
 * @param b
 * @return a∙b
 */
real dot_product(const Vector3& a,const Vector3& b);



/**
 * \~german Länge
 * \~english
 * \~
 * @param a
 * @return |a|
 */
real length(const Vector3& a);



/**
 * \~german
 * Vektor um die x Achse drehen
 *
 * \~english
 * Rotate vector on the x axis
 *
 * \~
 * @param p
 * @param a angle in rad
 * @return
 */
Vector3 rotate_x(const Vector3& p,real a);
Vector3 rotate_y(const Vector3& p,real a);
Vector3 rotate_z(const Vector3& p,real a);



//-------------------------
// spherical coordiantes
//-------------------------

/**
 * \~german
 * 3D-Polarkoordinaten (Kugelkoordinaten)
 * 
 * \~english
 * 3D Polar Coordinates (spherical coordinates)
 */
struct Polar3
{
	real r; ///< radial distance
	real th; ///< polar angle
	real ph; ///< azimuthal angle
};

/**
 * \~
 * Vector3 -> Polar3
 * 
 * @param v
 * @return 
 */
Polar3 polar3(const Vector3& v);

/**
 * \~
 * Polar3 -> Vector3
 * 
 * @param p
 * @return 
 */
Vector3 vector3(const Polar3& p);



//----------------
// functions
//----------------

/**
 * \~
 * @return a²-b²
 */
template<typename T>
constexpr T sub_sqr(T a,T b)
{
	return (a-b)*(a+b); // less rounding error
}


/**
 * \~
 * @return 1-x²
 */
template<typename T>
constexpr T _1_minus_sqr(T x)
{
	return sub_sqr(T(1),x);
}


/**
 * \~
 * @return √(1-x²)
 */
inline real sqrt_1_minus_sqr(real x)
{
	return std::sqrt(_1_minus_sqr(x));
}


/**
 * \~
 * @return x²
 */
template<typename T>
constexpr T sqr(T x)
{
	return x*x;
}



/**
 * \~
 * Signum
 * @return -1 or 0 or 1
 */
template <typename T>
int sgn(const T& x) {
	return (T(0) < x) - (x < T(0));
}




template<typename Iter,typename T>
T sample_variance(Iter it, Iter end,T mean)
{
	using I = typename std::iterator_traits<Iter>::value_type;
	static_assert(std::is_floating_point<I>::value,"floating point needed");
	static_assert(std::is_same<T,I>::value,"must be same type");
	
	const auto n = std::distance(it,end);		
	if (n<2)
		return T(NAN);
	
	I s = 0;
	for (;it!=end;++it) {		
		auto r = *it - mean;
		s += r*r;	
	}
	return s/(n-1);
}



template<typename T,size_t N>
T sample_variance(const std::array<T,N>& data,T mean)
{
	return sample_variance(data.begin(),data.end(),mean);
}



template<typename Iter,typename T = typename std::iterator_traits<Iter>::value_type>
T sample_standard_deviation(Iter begin, Iter end,T mean)
{	
	return std::sqrt(sample_variance(begin,end,mean));
}



template<typename T,size_t N>
T sample_standard_deviation(const std::array<T,N>& data,T mean)
{
	return sample_standard_deviation(data.begin(),data.end(),mean);
}



/**
 * \~german
 * Chi-Quadrat-Test (χ²-Test)
 *
 * \~english
 * χ² test
 *
 * \~
 * @param data
 * @param uncertainty
 * @param theory
 * @return
 */
template<typename T,size_t N>
T chi_squared_test(const std::array<T,N>& data,const std::array<T,N>& uncertainty,
				   const std::array<T,N>& theory,const size_t n=N)
{
	T chi = 0;
	auto di = data.begin();
	auto ui = uncertainty.begin();
	auto ti = theory.begin();
	auto end = data.begin() + n;
	for (;di!=end; ++di,++ui,++ti) {
		chi += sqr((*di-*ti)/ *ui);
	}
	return chi;
}



/**
 * \~german
 * Einen Wert x auf das Intervall [0,p)  normieren.
 * 
 * \~english
 * Normalise the value x to the interval [0,p).
 * 
 * \~
 * @param x value
 * @param p period
 * @return
 */
template<typename T>
T norm_period(T x,const T p)
{
	if (x < T(0)) {
		if (x < -p)
			x = std::fmod(x,p);
		x += p;
	}
	else {
		if (x >= p) {
			x -= p;
			if (x >= p)
				x = std::fmod(x,p);
		}
	}

	return x;
}



/**
 * \~german
 * Periodischer Abstand
 * 
 * Beispiel: 
 * Bei einer Periode von 1 ist der Abstand (absolute Differenz) zwischen 0,9 und 0,1 
 * gleich 0,2.
 * 
 * \~english
 * Example:
 * With a period of 1, the distance (absolute difference) between 0.9 and 0.1
 * equals 0.2.
 * 
 * \~
 * @param a
 * @param b
 * @param period
 * @return 
 */
double periodic_distance(double a,double b,double period);



/**
 * \~german
 * Wert auf das Intervall [0,360) normieren.
 *
 * \~english
 * Normalize value to the interval [0,360).
 *
 * \~
 * @param deg
 * @return
 */
real period_360(real deg);



/**
 * \~german
 * Wert auf das Intervall [0,2π) normieren.
 *
 * \~english
 * Normalize value to the interval [0,2π).
 *
 * \~
 * @param rad
 * @return
 */
real period_2pi(real rad);



/**
 * \~german
 * Wert auf das Interval [0,24) normieren.
 *
 * \~english
 * Normalize value to the interval [0,24).
 *
 * @param h
 * @return
 */
real period_24(real h);


/**
 * \~german a ungefähr gleich b?
 * Die ungefähre Gleichheit wird über das maximale Verhältnis bestimmt.
 *
 * \~english a approximate equal b?
 * Approximate equality is defined by the maximum ratio.
 *
 * \~
 * @param a
 * @param b
 * @param ratio
 * @return a≈b ?
 */
bool approximate(real a, real b,real ratio);



/**
 * \~german a ungefähr gleich b?
 *
 * Die ungefähre Gleichheit wird über die Differenz bestimmt.
 *
 * \~english a approximate equal b?
 *
 * Approximate equality is defined by the difference.
 *
 * \~
 * @param a
 * @param b
 * @param delta
 * @return a≈b ?
 */
bool approximate_delta(real a, real b,real delta);



/**
 * \~german a ungefähr gleich b?
 *
 * Wie approximate_delta(), aber unter Beachtung einer periodischen Funktion von a und b.
 *
 * \~english a approximate equal b?
 *
 * Like approximate_delta(), but taking into account a periodic function of a and b.
 *
 * \~
 * @param a
 * @param b
 * @param delta
 * @param period
 * @return a≈b ?
 */
bool approximate_delta_periodic(real a, real b,real delta,real period);



/**
 * \~german
 * Chi-Quadrat-Wahrscheinlichkeitsdichtefunktion
 *
 * \~english
 * chi squared probability density function
 *
 * \~
 * @param X χ²
 * @param n degrees of freedom
 * @return
 */
real chi_square_pdf(real X,int n);



/**
 * \~german
 * Chi-Quadrat-Verteilungsfunktion
 *
 * \~english
 * chi squared cumulative distribution function
 *
 * \~
 * @param X χ²
 * @param n degrees of freedom
 * @return
 */
real chi_square_cdf(real X,int n);



/**
  \~german
 * Intervallhalbierung
 *
 * Es ist erlaubt: x1>x2.
 * f()=0: Abbruch, f()>0: Richtung x2, f()<0: Richtung x1
 *
 * \~english
 * Bisection method
 * 
 * It is allowed: x1>x2.
 * f()=0: Cancelation, f()>0: direction x2, f()<0: direction x1
 * \~
 * @param x1
 * @param x2
 * @param f f(x1,p,x2) with p=current position
 *
 */
template<typename T,typename F>
void bisection(T x1,T x2,F f)
{
	static_assert(std::is_floating_point<T>::value,"floating point expected");

	auto p = (x1+x2)*T(0.5);
	while (p!=x1 && p!=x2) {
		auto r = f(x1,p,x2);
		if (r==T(0))
			break;
		(r>T(0) ? x1 : x2) = p;

		p = (x1+x2)*T(0.5);
	}	
}


/**
 * \~german
 * 1. Ableitung bestimmen
 *
 * Mathematisch genau bis Polynome 2. Grades.
 *
 * \~english
 * 1. derivative
 *
 * Mathematical exact until polynomials of 2nd degree.
 *
 * \~
 * @param x
 * @param f f(x)
 * @param h step
 * @return f'(x)
 */
template<typename T,typename F>
auto derive(T x,F f,T h) -> decltype(f(x))
{
	static_assert(std::is_floating_point<T>::value,"floating point type expected");

	return (f(x+h)-f(x-h))/(h+h);
}



/**
 * \~german
 * Gradient
 *
 * \~english
 * Gradient
 *
 * \~
 * @param f f(args)
 * @param[in][out] args arguments of f. Modified during call!
 * @param h step for derive()
 * @return gradient
 */
template<typename T,typename X,typename F>
X grad(F f,X& args,T h)
{
	X g;

	auto ai = std::begin(args);
	for(auto gi=std::begin(g); gi!=std::end(g); ++gi,++ai) {
		*gi = derive(*ai,[&ai,&f,&args](T x) {
			auto ax = *ai;
			*ai = x;
			auto y = f(args);
			*ai = ax;
			return y;
		},h);
	}

	return g;
}



template<typename X,typename T=typename X::value_type>
T vector_length(const X& x)
{
	T s=0;
	for (auto&& e : x) {
		s += sqr(e);
	}
	return std::sqrt(s);
}


template<typename X,typename T=typename X::value_type>
X normalize_vector(const X& x) {
	X y {x};

	T len = vector_length(x);
	for (auto&& e : y) {
		e /= len;
	}

	return y;
}





using Romberg_Cancelation_t = std::function<bool(int i,const std::vector<real>& Tk)>;

class RombergContinue
{
public:
	bool operator()(int,const std::vector<real>&) const {
		return false;
	}
};



/**
 * \~german
 * Integral von f(x) mit Romberg-Verfahren.
 *
 * \~english
 * Integral of f(x) using Romberg's method.
 *
 * \~
 * @param a
 * @param b
 * @param f f(x)
 * @param n 2^n pieces
 * @param cancel callback function to cancel the integration at some point
 * @return
 */
real integrate(real a,real b,std::function<real(real x)> f,int n,
			   Romberg_Cancelation_t cancel = RombergContinue());



/**
 * \~german
 * Dichtefunktion der Normalverteilung
 * 
 * \~english
 * Probability density function of the normal distribution
 * 
 * \~
 * @param m mu
 * @param s sigma
 * @param x
 * @return 
 */
real normal_pdf(real m,real s,real x);



/**
 * \~german
 * Verteilungsfunktion der Normalverteilung
 *
 * \~english
 * Cumulative distribution function for the normal distribution
 *
 * \~
 * @param m mu
 * @param s sigma
 * @param x
 * @return
 */
real normal_cdf(real m,real s,real x);





enum class ADTestType {
	Stephens,DAgostino
};



struct ADTestQuantiles {
	double q_50,q_25,q_15,q_10,q_5,q_1; // upper tail
};



const ADTestQuantiles ADTestStephensQuantiles  {  NAN,   NAN, 0.576, 0.656, 0.787, 1.092};
const ADTestQuantiles ADTestDAgostinoQuantiles {0.341, 0.470, 0.561, 0.631, 0.752, 1.035}; // isbn:0824774876, p. 123, Table 4.7



const ADTestQuantiles& test_quantiles(ADTestType type);



/**
 * \~german
 * Test auf Normalverteilung
 *
 * Anderson-Darling-Test wie gut eine Normalverteilung besteht.
 * Ohne Kenntnis der Parameter µ und σ der Normalverteilung.
 *
 * \~english
 * Anderson-Darling test for normality.
 * Without knowing µ and σ of the normal distribution.
 *
 * \~
 * @param sorted_samples Must be sorted ascending!
 * @return A²*
 */
double test_for_normality(const std::vector<double>& sorted_samples,
						  ADTestType transformation = ADTestType::DAgostino);




/**
 *\~german 
 * Konfidenzintervall für Binomialverteilungen
 * 
 * \~english
 * Binomial proportion Confidence interval 
 */
struct ConfidenceInterval
{
	double p; ///< success probability estimate
	double u; ///< uncertainty
};



/**
 * \~german
 * Agresti-Coull-Intervall für Binomialverteilungen.
 *
 * Für eine Konfidenz von 95%.
 *
 * \~english
 * Agresti-Coull interval for binomial distributions.
 *
 * For a 95% confidence level.
 * \~
 * @param k
 * @param n
 * @return
 */
ConfidenceInterval confidence_interval(int k,int n);



/**
 * \~german
 * p-Wert für die Größe A²* eines AD-Tests.
 * Die Genauigkeit ist begrenzt. Abweichungen im Bereich von 0.005 sind üblich.
 * (Nur für D'Agostino implementiert)
 *
 * \~english
 * p-value for the A²* value of an AD test.
 * The precision is limited. Deviations in the region of 0.005 are common.
 * (Only implemented for D'Agostino)
 *
 * \~
 * @param A A²*
 * @return
 */
double p_value_A(double A,ADTestType type = ADTestType::DAgostino);



}//aether


#endif // MATHEMATICS_H
