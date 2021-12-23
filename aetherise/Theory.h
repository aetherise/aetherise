/**
 * \file
 *
 * \~german
 * @brief Theorieschnittstelle und Implementierungen
 *
 * \~english
 * @brief Theory interface and implementations
 *
 */

#ifndef THEORY_H
#define THEORY_H

#include "mathematics.h"
#include "utils.h"
#include "physics.h"

#include <vector>

namespace aether {



/**
 * \~german Lichtweg
 *
 * Eine Strecke mit Brechungsindex des Materials.
 * Der Vektor p gibt die Richtung und die Länge der Strecke an.
 *
 * \~english
 * A line segment with index of refraction of the material.
 * The vector p defines the direction and the length of the line segment.
 */
struct Lightpath
{
	Vector3 p; ///< direction and length [m]
	double n; ///< index of refraction
};



/**
 * \~german
 * Ein Interferometer mit einer Lichtquelle und zwei Lichtwegen,
 * wie bei einem Michelson-Interferometer.
 *
 * \~english
 * An interferometer with one light source and two light paths,
 * like a Michelson interferometer.
 */
struct Interferometer
{
	double wave_length; // [m]
	std::vector<Lightpath> light_path_1;
	std::vector<Lightpath> light_path_2;
};



/**
 * \~german
 * Theorie (Schnittstelle)
 *
 * \~english
 * Theory (interface)
 */
class Theory
{
public:
	virtual ~Theory() = default;

	/**
	 * \~german 
	 * Phasenlaufzeit
	 * 
	 * Reisezeit eines Photons für die gegebene Strecke.
	 *
	 * \~english
	 * Travel time of a photon for the given line segment.
	 *
	 * \~
	 * @param v velocity vector of movement relative to the aether
	 * @param light_path
	 * @return duration in [s]
	 */
	virtual real phase_delay(const Vector3& v,const Lightpath& light_path) const = 0;
};



/**
 * \~german
 * Klassische Äthertheorie
 *
 * Keine Lorentz-Kontraktion und kein anisotroper Brechungsindex bei Gasen.
 *
 * \~english
 * Classic aether theory
 *
 * No Lorentz contraction and no anisotropic refractive index of gases.
 */
class Classic : public Theory
{
public:
	real phase_delay(const Vector3& v,const Lightpath& light_path) const override;
};





/**
 * \~german
 * Äthertheorie von Lorentz
 *
 * \~english
 * Lorentz aether theory
 */
class Aether : public Theory
{
public:
	real phase_delay(const Vector3& v,const Lightpath& light_path) const override;
};



/**
 * \~german
 * Spezielle Relativitätstheorie
 *
 * \~english
 * Special Theory of Relativity
 */
class Relativity : public Theory
{
public:
	real phase_delay(const Vector3& v,const Lightpath& light_path) const override;
};




/**
 * \~german
 * Phasenlaufzeit eines Photons für einen Lichtweg bestehend aus mehreren Strecken berechnen.
 *
 * \~english
 * Phase delay of a photon on a light path consisting of multiple line segments.
 *
 * \~
 * @param v
 * @param light_path
 * @param theory
 * @return duration in [s]
 */
real phase_delay(const Theory& theory,const Vector3& v,const std::vector<Lightpath>& light_path);



/**
 * \~german Streifenverschiebung
 *
 * Streifenverschiebungen von Millers Interferometer auf dem Mount Wilson.
 * Berechnet wird die Verschiebung in Wellenlängen. Eine Verschiebung von 0,5 bedeutet,
 * daß helle und dunkle Streifen die Plätze tauschen.
 *
 * \~english
 * Fringe displacements of Miller's interferometer on the Mount Wilson.
 * Calculates the displacement in wave length. A displacement of 0.5 means
 * bright and dark fringes exchange positions.
 *
 * \~
 * @param theory
 * @param params velocity vector of movement in the aether
 * @param latitude in rad
 * @param n index of refraction
 * @param L effective arm length
 * @param sidereal_time in rad
 * @param invert change sign of values?
 * @return 17 samples of a complete rotation
 */
std::array<double,17>
fringe_displacements(const Theory& theory,const TheoryParameters& params,double latitude, double n,double L,
					 double sidereal_time,bool invert);



/**
 * \~german
 * Bestimmtes neues Theorieobjekt erzeugen.
 *
 * \~english
 * Create a new specific theory object.
 *
 * \~
 * @param options_theory
 * @return
 */
std::unique_ptr<Theory> create_theory(Options::Theory options_theory);




}//aether

#endif // theory_H
