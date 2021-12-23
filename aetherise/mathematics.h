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

#include "stdx.h"

#include <cmath>
#include <numeric>
#include <stdexcept>
#include <vector>
#include <functional>
#include <limits>
#include <complex>
#include <map>

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


constexpr double rad(double deg) {
	return AETHER_PI/180.0 * deg;
}

constexpr double deg(double rad) {
	return 180.0/AETHER_PI * rad;
}

constexpr double rad_to_h(double rad) {
	return 12.0/AETHER_PI*rad;
}

constexpr double h_to_rad(double h) {
	return AETHER_PI/12.0*h;
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
template<typename Iter,typename F>
auto mean_value(Iter it,Iter end,F get) -> decltype(get(*it))
{	
	using T = decltype(get(*it));	
	static_assert(std::is_floating_point<T>::value,"floating point needed");

	auto n = std::distance(it,end);
	if (n <= 0)
		return T(NAN);

	T mean = get(*it);
	for (++it; it != end; ++it)
		mean += get(*it);
	return mean/T(n); // T/unsigned !
}



template<typename Iter,typename T=typename std::iterator_traits<Iter>::value_type>
T mean_value(Iter it,Iter end)
{
	return mean_value(it,end,[](const T& value){
		return value;
	});
}



template<typename T,size_t N>
T mean_value(const std::array<T,N>& a)
{
	return mean_value(a.begin(),a.end());
}


/**
 * \~german Median
 * Der Wert des Elementes in der Mitte.
 * Bei gerader Anzahl an Werten wird der Mittelwert der beiden Werte 
 * in der Mitte berechnet.
 * 
 * \~english
 * The value of the element in the middle.
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
template<typename T>
T sqrt_1_minus_sqr(T x)
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



/**
 * \~german 
 * Stichprobenvarianz (erwartungstreu)
 * 
 * 
 * \~english
 * Sample variance (corrected)
 * 
 * \~
 * @param mean estimate of expected value
 */
template<typename Iter,typename T,typename F>
auto sample_variance(Iter it, Iter end,T mean,F get) -> decltype (get(*it))
{
	using G = decltype (get(*it));
	static_assert(std::is_same<T,G>::value,"must be same type");
	static_assert(std::is_floating_point<T>::value,"floating point needed");
	
	
	const auto n = std::distance(it,end);		
	if (n<2)
		return T(NAN);
	
	T s = 0;
	for (;it!=end;++it) {				
		s += sqr(get(*it) - mean);	
	}
	return s/T(n-1);
}



template<typename Iter,typename T>
T sample_variance(Iter it, Iter end,T mean)
{		
	return sample_variance(it,end,mean,[](const T& sample){
		return sample;
	});
}



/**
 * \~german
 * Stichprobenstandardabweichung
 * 
 * \~english
 * Experimental standard deviation (GUM)
 * 
 * \~
 * @param mean estimate of expected value
 */
template<typename Iter,typename T>
T sample_standard_deviation(Iter begin, Iter end,T mean)
{	
	return std::sqrt(sample_variance(begin,end,mean));
}



/**
 * \~german
 * Residuen
 * 
 * \~english
 * 
 * \~
 * @param data
 * @param uncertainty
 * @param theory
 * @param n
 * @return 
 */
template<typename T,size_t N>
std::array<T,N> residuals(const std::array<T,N>& data,const std::array<T,N>& uncertainty,
						  const std::array<T,N>& theory,const size_t n=N)
{
	std::array<T,N> r;
	auto di = data.begin();
	auto ui = uncertainty.begin();
	auto ti = theory.begin();
	auto ri = r.begin();
	auto end = data.begin() + n;
	for (;di!=end; ++di,++ui,++ti,++ri) {
		*ri = (*di-*ti)/ *ui;
	}
	return r;
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
 * Normalize the value x to the interval [0,p).
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
 * \~german Phasendifferenz
 * 
 * \~english
 * 
 * \~
 * @param p0 reference
 * @param p
 * @param period
 * @return 
 */
double phase_difference(double p0,double p,double period);



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
double period_360(double deg);



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
double period_2pi(double rad);



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
double period_24(double h);


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
bool approximate(double a, double b,double ratio);



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
bool approximate_delta(double a, double b,double delta);



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
bool approximate_delta_periodic(double a, double b,double delta,double period);



/**
 * \~german
 * Chi-Quadrat-Wahrscheinlichkeitsdichtefunktion
 *
 * \~english
 * chi squared probability density function
 *
 * \~
 * \warning if type of real is double, precision loss at n>171
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
 * \warning if type of real is double, precision loss at n>171
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
auto differentiate(T x,F f,T h) -> decltype(f(x))
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
 * @param h step for differentiate()
 * @return gradient
 */
template<typename T,typename X,typename F>
X grad(F f,X& args,T h)
{
	X g;

	auto ai = std::begin(args);
	for(auto gi=std::begin(g); gi!=std::end(g); ++gi,++ai) {
		*gi = differentiate(*ai,[&ai,&f,&args](T x) {
			auto ax = *ai;
			*ai = x;
			auto y = f(args);
			*ai = ax;
			return y;
		},h);
	}

	return g;
}




using Romberg_Cancelation_t = std::function<bool(int i,const std::vector<double>& Tk)>;

class RombergContinue
{
public:
	bool operator()(int,const std::vector<double>&) const {
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
double integrate(double a,double b,std::function<double(double x)> f,int n,
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



/**
 * \~german
 * Testquantile der Größe A²*. 
 * Ist ein Wert größer als ein Quantil, ist die Normalverteilung abgelehnt.
 * 
 * \~english
 * Test quantiles of A²*.
 * If a value is larger than a quantile, the normality is rejected.
 * 
 * \~
 * @param type
 * @return 
 */
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
 * @param sorted_samples Must be sorted ascending! n>=8.
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




/**
 * \~german
 * Ausreißer erkennen mit Chauvenets Kriterium.
 * Liefert den schlimmsten Ausreißer zurück.
 * 
 * \~english
 * Detect outliers with Chauvenet's criterion.
 * Finds the worst outlier.
 * 
 * \~
 * TODO not well tested
 * @param begin
 * @param end
 * @return end if no outlier detected
 */
template<typename Iter>
Iter detect_outlier(Iter begin, Iter end)
{		
	using T = typename std::iterator_traits<Iter>::value_type;
	static_assert (std::is_floating_point<T>::value,"Floating point expected");

	const auto n = std::distance(begin,end);		
	if (n<4)
		return end;
	
	T mean = mean_value(begin,end);		
	T s = sample_standard_deviation(begin,end,mean);
	
	const auto P = T(1)-T(0.5)/n; // Chauvenet's criterion
	auto max_it = end;
	T max_p = 0;
	for (Iter it = begin;it!=end;++it) {		
		auto k = std::abs(*it-mean)/s;		
		auto p = std::erf(k*T(1./std::sqrt(2.)));
		if (p > P) {
			if (p > max_p) {
				max_p = p;
				max_it = it;
			}
		}
	}		
		
	return max_it;
}



/**
 * \~german
 * Diskrete Fourier Transformation (DFT) - Goertzel Algorithmus
 * 
 * \~english
 * Discrete Fourier Transform (DFT) - Goertzel Algorithm
 *  
 * \~
 * @param k harmonic
 * @param it iterator
 * @param end iterator
 * @return Amplitude and phase of a cosine
 */
template <typename Iter>
std::complex<double> DFT_analyze(typename std::iterator_traits<Iter>::difference_type k,Iter it, Iter end)
{	
	// Gerald Goertzel, AN ALGORITHM FOR THE EVALUATION OF FINITE TRIGONOMETRIC SERIES
	// The American Mathematical Monthly, Vol. 65, No. 1 (Jan., 1958), pp. 34-35
		
	using K = typename std::iterator_traits<Iter>::difference_type;
	
	auto n = std::distance(it,end);		
	if (n<=0)
		return {};
	k = clamp(k,K(0),n-1);
			
	double w = (AETHER_2PI/n)*k;
	double cosw = std::cos(w)*2;	
	double sinw = std::sin(w);
	
	double U1 = 0;
	double U2 = 0;
			
	for(; it != end; ++it) {
		double Uk = *it + cosw*U1 - U2;		
		U2 = U1;
		U1 = Uk;				
	}
	// one iteration more, to get the correct phase
	double Uk = cosw*U1 - U2;		
	U2 = U1;
	U1 = Uk;
	
	double C = U1 - U2*0.5*cosw; 
	double S = U2*sinw;
	C = C/n*2;
	S = S/n*2; 
	
	return {C,S};
}



/**
 *\~german
 * Diskrete Fourier Transformation (DFT) - Goertzel Algorithmus
 * 
 * \~english
 * Discrete Fourier Transform (DFT) - Goertzel Algorithm
 */
class DFTGoertzel
{
	struct State {		
		double cosw,sinw;
		double U1,U2;
		
		std::complex<double> result(size_t n) const;		
	};
	
	std::vector<double> frequencies;	
	std::vector<State> states;
	size_t _n;
	
public:
	/**
	 * \~german
	 * 
	 * \~english
	 * 
	 * \~
	 * @param frequencies 
	 * @param sample_rate Number of samples for frequency 1
	 */
	DFTGoertzel(const std::vector<double>& frequencies,int sample_rate);
	  
	
	/**
	 *
	 */
	template<typename Iter>	
	void analyze(Iter it,Iter end) {
		auto n = std::distance(it,end);		
		if (n<=0)
			return;

		_n += n;
		
		for(; it != end; ++it) { // possible long array loop outside, better performance?
			for (auto& state : states) {		
				double Uk = *it + state.cosw*state.U1 - state.U2;		
				state.U2 = state.U1;
				state.U1 = Uk;		
			}		
		}
	}
	
	
	/**
	 *
	 * @return 
	 */
	std::map<double,std::complex<double>> result() const;
	
};




template<typename T,size_t N>
void set_sine(double A, double phi,double c,int sample_rate, std::array<T,N>& samples)
{	
	for (size_t i=0;i<N;i++) {
		samples[i] = A*std::sin(i*AETHER_2PI/sample_rate + phi) + c;		
	}
}



template<typename T,size_t N>
void add_sine(double A, double phi,double c,int sample_rate, std::array<T,N>& samples)
{	
	for (size_t i=0;i<N;i++) {
		samples[i] += A*std::sin(i*AETHER_2PI/sample_rate + phi) + c;		
	}
}



/**
 * \~german
 * Schätzwert mit beigeordneter Standardunsicherheit
 * 
 * \~english
 * Estimate with standard uncertainty
 * 
 */
template<typename T>
struct Estimate
{
	T m,u;
};



/**
 * \~german
 * Erwartungswert schätzen und Standardunsicherheit bestimmen.
 * Die Werte müssen normal verteilt sein.
 * 
 * \~english
 * Estimate the expected value and evaluate the standard uncertainty.
 * The values must be normal distributed.
 * 
 * \~
 * @param begin 
 * @param end
 * @param get Getter for the value of an element
 * @return Estimate with uncertainty
 */
template<typename Iter,typename F>
auto estimate_expected_value(Iter begin,Iter end,F get) -> Estimate<decltype(get(*begin))>
{			
	using T = decltype(get(*begin));
	const auto n = std::distance(begin,end);
	
	Estimate<T> e;			
	e.m = mean_value(begin,end,get);
	e.u = sample_variance(begin,end,e.m,get);			
	e.u = std::sqrt(e.u/n); // standard uncertainty
		
	return e;
}




/**
 * \~german
 * Erwartungswert schätzen und Standardunsicherheit bestimmen.
 * Die Werte müssen normal verteilt sein.
 * 
 * \~english
 * Estimate the expected value and evaluate the standard uncertainty.
 * The values must be normal distributed.
 * 
 * \~
 * @param begin 
 * @param end  
 * @return Estimate with uncertainty
 */
template<typename Iter,typename T=typename std::iterator_traits<Iter>::value_type>
Estimate<T> estimate_expected_value(Iter begin,Iter end)
{	
	return estimate_expected_value(begin,end,[](const T& sample){
		return sample;
	});
}




/**
 * \~german
 * Realteil mit Unsicherheit
 * 
 * \~english 
 * Real part with uncertainty
 * 
 * \~
 * @param z 
 * @return
 */
template<typename T>
Estimate<T> Real(const Estimate<std::complex<T>>& z)
{
	return {z.m.real(),z.u.real()};
}



/**
 * \~german
 * Imaginärteil mit Unsicherheit
 * 
 * \~english 
 * Imaginery part with uncertainty
 * 
 * \~
 * @param z 
 * @return
 */
template<typename T>
Estimate<T> Imag(const Estimate<std::complex<T>>& z)
{
	return {z.m.imag(),z.u.imag()};
}



/**
 *\~
 * @return a+b
 */
template<typename T>
Estimate<T> propagate_add(const Estimate<T>& a, const Estimate<T>& b)
{
	return {a.m+b.m,std::sqrt(sqr(a.u)+sqr(b.u))};
}



/**
 * \~
 * @return a-b
 */
template<typename T>
Estimate<T> propagate_sub(const Estimate<T>& a, const Estimate<T>& b)
{
	return {a.m-b.m,std::sqrt(sqr(a.u)+sqr(b.u))};
}



/**
 *\~
 * @return a+b
 */
template<typename T>
Estimate<std::complex<T>> propagate_add(const Estimate<std::complex<T>>& a, const Estimate<std::complex<T>>& b)
{	
	auto uR = std::sqrt(sqr(a.u.real())+sqr(b.u.real()));
	auto uI = std::sqrt(sqr(a.u.imag())+sqr(b.u.imag()));
	return {a.m+b.m,{uR,uI}};
}



/**
 * \~
 * @return a-b
 */
template<typename T>
Estimate<std::complex<T>> propagate_sub(const Estimate<std::complex<T>>& a, const Estimate<std::complex<T>>& b)
{
	auto uR = std::sqrt(sqr(a.u.real())+sqr(b.u.real()));
	auto uI = std::sqrt(sqr(a.u.imag())+sqr(b.u.imag()));
	return {a.m-b.m,{uR,uI}};
}





/**
 * \~german Stichprobenkovarianz
 * 
 * Kovarianz der Komponenten von komplexen Zahlen.
 * 
 * \~english
 * Covariance of the components of complex numbers.
 * 
 * \~
 * @param e expected value
 * @param it begin
 * @param end end
 * @return covariance
 */
template<typename T, typename Iter>
T sample_covariance(const std::complex<T>& e,Iter it,Iter end)
{
	// TODO test
	auto n = std::distance(it,end);
	if (n<=0)
		return NAN;
	
	T s = 0;
	for (; it!=end; ++it) {				
		s += (it->real() - e.real())*(it->imag() - e.imag());	
	}
	s /= T(n-1);
		
	return s;
}





/**
 * \~german
 * Fehlerfortpflanzung durch Funktion f()
 * 
 * \~english
 * Error propagation through function f()
 * 
 * \~
 * @param e Argument for f()
 * @param f Function y=f(e)
 * @param h step parameter for numeric derivative
 * @return y with uncertainties
 */
template<typename T,typename F>
Estimate<T> propagate_uncertainties(const Estimate<T>& e,F f,T h)
{
	auto dfde = differentiate(e.m,[&f](const T& e){
		return f(e);
	},h);
	return {f(e.m),e.u*dfde};
}



/**
 * \~german
 * Fehlerfortpflanzung durch Funktion f()
 * 
 * \~english
 * Error propagation through function f()
 * 
 * \~
 * @param e1 First argument for f()
 * @param e2 Second argument for f()
 * @param f Function y=f(e1,e2)
 * @param h step parameter for numeric derivative
 * @param cov estimated covariance = sample_cov(e1,e2)/n
 * @return y with uncertainties
 */
template<typename T,typename F>
Estimate<T> propagate_uncertainties(const Estimate<T>& e1,const Estimate<T>& e2,F f,T h,T cov=0)
{	
	// partial derivatives
	auto dfde1 = differentiate(e1.m,[&e2,&f](const T& e1){
		return f(e1,e2.m);
	},h);
	auto dfde2 = differentiate(e2.m,[&e1,&f](const T& e2){
		return f(e1.m,e2);
	},h);
	
	auto s = sqr(e1.u*dfde1)+sqr(e2.u*dfde2) + 2*dfde1*dfde2*cov;
	if (s<0) // should not happen, but can occur because of rounding errors
		s = 0;
	auto u = std::sqrt(s);
	return {f(e1.m,e2.m),u};
}




}//aether


#endif // MATHEMATICS_H
