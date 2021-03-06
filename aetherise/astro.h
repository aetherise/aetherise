/**
 * \file
 *
 * \~german
 * @brief Astronomie
 *
 * \~english
 * @brief astronomy
 *
 */


#ifndef ASTRO_H
#define ASTRO_H

#include "mathematics.h"
#include "utils.h"


namespace aether {

constexpr const double Earth_average_orbital_speed = 29780; // (m/s)

constexpr const double Epoch_J2000 = 2451545.0;




/**
 * \~german Julianisches Jahrhundert
 * 
 * \~english
 * 
 * \~
 * @param JD
 * @param epoch
 * @return 
 */
constexpr double julian_century(double JD,double epoch)
{
	return (JD-epoch)/36525;
}



/**
 * \~german
 * Schiefe der Ekliptik.
 *
 * \~english
 *
 * \~
 * @param JD julian day
 * @return
 */
double obliquity_of_the_ecliptic(double JD);



/**
 * \~german
 * Horizontales Koordinatensystem
 *
 * \~english
 * Horizontal coordinate system
 */
struct Horizontal
{
	real a; ///< azimuth (from south, π/2=west) in (rad)
	real h; ///< height in (rad)
};


/**
 * \~german
 * Äquatorial ruhendes Koordinatensystem
 *
 * \~english
 * Equatorial coordinate system at rest
 */
struct Equatorial_resting
{
	real t; ///< hour angle in (rad)
	real d; ///< declination in (rad)
};



/**
 * \~german
 * Äquatorial rotierendes Koordinatensystem
 *
 * \~english
 * equatorial coordinate system
 */
struct Equatorial
{
	real ra; ///< right ascension in (rad)
	real de; ///< declination in (rad)
};



/**
 * \~german
 * Ekliptikales Koordinatensystem
 * (Geozentrisch)
 *
 * \~english
 * Ecliptic coordinate system
 * (geocentric)
 */
struct Ecliptic
{
	real l; ///< longitude in (rad)
	real b; ///< latitude in (rad)
};


/**
 * \~german
 * Galaktisches Koordinatensystem
 *
 * \~english
 * Galactic coordinate system
 *
 */
struct Galactic
{
	real l; ///< longitude in (rad)
	real b; ///< latitude in (rad)
};



/**
 * \~german
 * Äquatorial -> Ekliptikal
 *
 * \~english
 * equatorial -> ecliptic
 *
 * \~
 * @param eq
 * @param obliquity of the ecliptic
 * @return
 */
Ecliptic ecliptic(const Equatorial& eq, const real e);



/**
 * \~german
 * Ekliptikal -> Äquatorial
 *
 * \~english
 * ecliptic -> equatorial
 *
 * \~
 * @param ec
 * @param obliquity of the ecliptic
 * @return
 */
Equatorial equatorial(const Ecliptic& ec, const real e);



/**
 * \~german
 * Äquatorial ruhend -> rotierend
 *
 * \~english
 * equatorial resting -> rotating
 *
 * \~
 * @param er
 * @param sidereal_time in (rad)
 * @return
 */
Equatorial rotating(const Equatorial_resting& er,real sidereal_time);



/**
 * \~german
 * Äquatorial rotierend -> ruhend
 *
 * \~english
 * equatorial rotating -> resting
 *
 * \~
 * @param eq
 * @param sidereal_time in (rad)
 * @return
 */
Equatorial_resting resting(const Equatorial& eq,real sidereal_time);


/**
 * \~german
 * Horizontal -> Äquatorial ruhend
 *
 * \~english
 * horizontal -> equatorial resting
 *
 * \~
 * @param hz
 * @param latitude in (rad)
 * @return
 */
Equatorial_resting equatorial(const Horizontal& hz,real latitude);


/**
 * \~german
 * Äquatorial ruhend -> Horizontal
 *
 * \~english
 * equatorial resting -> horizontal
 *
 * \~
 * @param er
 * @param latitude in (rad)
 * @return
 */
Horizontal horizontal(const Equatorial_resting& er,real latitude);


/**
 * \~german
 *  Galaktisch -> Äquatorial
 * 
 * \~english
 *  galactic -> equatorial
 * \~
 * @param gc
 * @return
 */
Equatorial equatorial(const Galactic& gc);


/**
 * \~german
 * Äquatorial -> Galaktisch
 * 
 * \~english
 * equatorial -> galactic
 * 
 * \~
 * @param ae
 * @return
 */
Galactic galactic(const Equatorial& eq);



/**
 * \~german
 * Südazimut -> Nordazimut
 *
 * \~english
 * south azimuth -> north azimuth
 *
 * \~
 * @param a south azimuth in (rad)
 * @return north azimuth
 */
constexpr real north_azimuth(real a)
{
	return a + AETHER_PI_L;
}





/**
 * \~german
 * Sternzeit
 *
 * \~english
 *
 * \~
 * @param cal date and time in UT (Universal Time)
 * @param longitude (rad) (negativ = west)
 * @return angle (rad)
 */
double sidereal_time(const Calendar& cal,double longitude);




/**
 * \~german
 * Koordinaten der Sonne (geozentrisch) zum gegebenen Zeitpunkt.
 *
 * \~english
 * Coordinates of the Sun (geocentric) for the given point in time.
 *
 * \~
 * @param JD julian day (UT)
 * @return
 */
Ecliptic sun_coordinates(double JD);



/**
 * \~german
 * Erdapex
 *
 * Die Richtung der Bewegung der Erde, auf ihrer Bahn um die Sonne,
 * zum angegebenen Zeitpunkt.
 *
 * \~english
 * The direction of movement of the earth, on its orbit around the sun,
 * for the given point in time.
 *
 * \~
 * @param JD julian day (UT)
 * @return
 */
Equatorial earth_apex(double JD);



}//aether

#endif // ASTRO_H
