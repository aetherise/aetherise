/**
 * \file
 *
 * \~german
 * @brief Minimierer
 *
 * \~english
 * @brief Minimizer
 *
 */

#ifndef MINI_H
#define MINI_H

#include "mathematics.h"

#include <ostream>
#include <vector>
#include <type_traits>


namespace aether {


/**
 * \~german Fortschritt nicht schätzen
 *     
 * \~english
 * 
 */
class DontEstimateProgress
{
public:
	void operator() (int ,double ,double ) {}
};



/**
 * \~german Lokales Minimum
 * 
 * Ergebnis einer Minimierung.
 * Wenn i gleich der maximalen Anzahl an Iterationen ist,
 * ist die Minimierung fehlgeschlagen.
 * 
 * \~english
 * Result of a minimization. The minimization failed, if i is
 * equal to the maximum number of iterations.
 */
template<typename X,typename Y>
struct LocalMinimum
{
	X x; ///< Parameters at minimum
	Y y; ///< f(x) at minimum
	int i; ///< Number of iterations needed
};



/**
 * \~german
 * Lokal minimieren (a)
 *
 * Findet ein lokales Minimum mit einem Gradientenverfahren.
 * Für Funktionen f:Rn->R1.
 *
 * \~english
 * Finds the local minimum using a gradient descent.
 * For functions f:Rn->R1.
 *
 * \~
 * \warning Experimental
 *
 * @param x0 starting point
 * @param f f(x)
 * @param h step for derive()
 * @param n maximum number of iterations
 * @param eps relative precision
 * @param progress p(i,y,y1)
 * @return
 */
template<typename T,typename X,typename F,typename ProgressF=DontEstimateProgress>
LocalMinimum<X,T>
minimize_locally_a(const X& x0,F f,T h,int n,T eps,ProgressF progress=ProgressF())
{
	static_assert(std::is_floating_point<T>::value,"floating point type expected");

	auto x = x0;
	auto y = f(x0);
	auto x1 = x;
	auto y1 = y;

	auto xmin = x;
	auto ymin = y;

	int i=0;
	while(i<n) {
		auto g = grad(f,x,h);
		
		// find x1
		auto gi = std::begin(g);
		auto x1i = std::begin(x1);
		for (auto xi=std::begin(x); xi!=std::end(x); ++xi,++gi,++x1i) {
			T gi0 = std::abs(*gi)<T(1) ? *gi : T(1)/ *gi;
			T c = h;
			T xi0 = *xi;
			T xc, yc;
			T xc1 = *xi;
			T yc1 = y;			
			do {
				xc = xc1;
				yc = yc1;
				auto dx = c*gi0;
				*xi = xi0 - dx;
				xc1 = *xi;
				yc1 = f(x);								
				c += c;				
			}
			while (yc1<=yc && c!=INFINITY);

			*x1i = xc;
			*xi = xi0;
		}

		y1 = f(x1);

		if (y1<ymin) {
			xmin = x1;
			ymin = y1;
		}

		auto dy = std::abs(y-y1);
		if (dy<=std::abs(y)*eps || (dy<eps && std::min(y,y1)==T(0)))
			break;	

		i++;
		progress(i,y,y1);
		x = x1;
		y = y1;
	}

	progress(n,y,y1);

	return {xmin,ymin,i};
}




/**
 * \~german
 * Lokal minimieren (b)
 *
 * Findet ein lokales Minimum mit einem Gradientenverfahren.
 * Für Funktionen f:Rn->R1.
 *
 * \~english
 * Finds the local minimum using a gradient descent.
 * For functions f:Rn->R1.
 *
 * \~
 * \warning Experimental
 *
 * @param x0 starting point
 * @param f f(x)
 * @param h step for derive()
 * @param n maximum number of iterations
 * @param eps relative precision
 * @param progress p(i,y,y1)
 * @return
 */
template<typename T,typename X,typename F,typename ProgressF=DontEstimateProgress>
LocalMinimum<X,T>
minimize_locally_b(const X& x0,F f,T h,int n,T eps,ProgressF progress=ProgressF())
{
	static_assert(std::is_floating_point<T>::value,"floating point expected");

	auto x = x0;
	auto y = f(x0);
	auto x1 = x;
	auto y1 = y;

	auto xmin = x;
	auto ymin = y;

	int i=0;
	while(i<n) {
		auto c = h;
		auto xc = x;
		auto yc = y;		
		T mdx;
		do {			
			x1 = xc;
			y1 = yc;			

			mdx = T(0);
			auto g = grad(f,xc,h);			
			auto gi = std::begin(g);
			for (auto xi=std::begin(xc); xi!=std::end(xc); ++xi,++gi) {				
				T gi0 = std::abs(*gi)<T(1) ? *gi : T(1)/ *gi;
				auto dx = gi0 * c;
				*xi = *xi - dx;
				mdx = std::max(mdx,std::abs(dx));
			}

			yc = f(xc);			
			c += c;
		}
		while((yc<=y || mdx<h) && c!=INFINITY);

		if (y1<ymin) {
			xmin = x1;
			ymin = y1;
		}

		auto dy = std::abs(y-y1);
		if (dy<=std::abs(y)*eps || (dy<eps && std::min(y,y1)==T(0)))
			break;

		i++;
		progress(i,y,y1);
		x = x1;
		y = y1;
	}

	progress(n,y,y1);

	return {xmin,ymin,i};
}



/**
 * \~german
 * Ergebnis einer Minimierung
 *
 * \~english
 * Result of a minimization
 */
struct MinimizerResult
{
	std::vector<double> x; ///< Parameters
	std::vector<double> u; ///< Uncertainties
	double y; ///< Function value at minimum
	bool valid; ///< Is result valid?
};



#ifdef AETHER_MINUIT

/**
 * \~german
 * Stellt alle Informationen bereit welche ein externer Minimierer benötigt.
 *
 * \~english
 * Provides all information which is needed by an external minimizer.
 */
class MinimizeContext
{
public:
	virtual ~MinimizeContext()=default;

	/**
	 * \~german
	 * Die Funktion die minimiert werden soll.
	 * 
	 * \~english
	 * The function to minimize.
	 * 
	 * \~
	 * @param params parameters
	 * @return 
	 */
	virtual double operator() (const std::vector<double>& params) = 0;

	
	
	/**
	 * \~german
	 * Über ∆χ² wird das Konfidenzniveau für die Unsicherheiten der Parameter am Minimum festgelegt.
	 * Ein Wert von 1 bedeutet die 1-Sigma Standardunsicherheit.
	 * 
	 * \~english
	 * Sets the confidence level of the uncertainties of the parameters at the minimum via the ∆χ².
	 * A value of 1 means 1-Sigma standard uncertainty.
	 * 
	 * \~
	 * @return ∆χ²
	 */
	virtual double delta_chi_squared() const = 0;
		
};



/**
 * \~german
 * Lokal minimieren (mit Minuit2)
 *
 * Nutzt ROOT::Minuit2 zur Minimierung.
 * Für Funktionen f:Rn->R1.
 *
 * \~english
 * Uses ROOT::Minuit2 to minimize.
 * For functions f:Rn->R1.
 *
 * \~
 * @param x0 start parameters
 * @param context information
 * @param os text output will be redirected here if possible
 * @return
 */
MinimizerResult
minimize_locally_Minuit2(const std::vector<double>& x0, MinimizeContext& context,std::ostream& os);

#endif // Minuit2



}//aether

#endif // MINI_H
