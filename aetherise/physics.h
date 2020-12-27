/**
 * \file
 *
 * \~german
 * @brief Physik
 *
 * \~english
 * @brief Physics
 *
 */


#ifndef PHYSICS_H
#define PHYSICS_H



namespace aether {


const double Lightspeed = 299792458.0; // (m/s)

const double Gas_Constant = 8.31446261815324; // (kg m²)/(s² mol K)


/**
 *
 * @param t °C
 * @return K
 */
constexpr double in_Kelvin(double t)
{
	return t + 273.15;
}

/**
 *
 * @param T K
 * @return °C
 */
constexpr double in_Celsius(double T)
{
	return T - 273.15;
}



/**
 * \~german
 * Brechungsindex von Luft
 *
 * Verwendet die Formeln von Philip E. Ciddor.
 *
 * @param p Druck in Pascal
 * @param T Temperatur in Kelvin
 * @param lambda Wellenlänge in µm
 * @param h relative Luftfeuchtigkeit [0,1]
 * @param xc CO2 Konzentration in ppm
 * @return
 *
 * \~english
 * Refractive index of air
 *
 * Based on Ciddor Equation.
 *
 * @param p pressure in Pascal
 * @param T temperature in Kelvin
 * @param lambda wave length in µm
 * @param h relative humidity [0,1]
 * @param xc CO2 content in ppm
 * @return
 */
double refractive_index_of_air(double p,double T,double lambda,double h,double xc);



/**
 * \~german
 * Barometrische Höhenformel
 *
 * zur Berechnung des Luftdrucks in der Höhe h.
 * Gilt nährungsweise bis in eine Höhe von ungefähr 11 km.
 * Die Parameter p0 und T0 sind auf die Werte der internationalen
 * Standardatmosphäre voreingestellt.
 *
 * \~english
 * Calculates the atmospheric pressure at the given height h.
 *
 * Approximately valid only to a height of 11 km.
 * Uses the international standard atmosphere
 * as default for the parameters p0 and T0.
 *
 * \~
 * @param h height in m
 * @param p0 in Pa
 * @param T0 in K
 * @return air pressure in Pa
 */
double barometric_formula(double h,double p0=101325,double T0=288.15);



}//


#endif
