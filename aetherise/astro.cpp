
#include "astro.h"

namespace aether {



Equatorial rotating(const Equatorial_resting& er,real sidereal_time)
{
	Equatorial eq;
	eq.ra = sidereal_time - er.t;
	eq.de = er.d;
	return eq;
}


Equatorial_resting resting(const Equatorial& eq,real sidereal_time)
{
	Equatorial_resting er;
	er.t = sidereal_time - eq.ra;
	er.d = eq.de;
	return er;
}


Equatorial_resting equatorial(const Horizontal& hz,real latitude)
{
	using namespace std;
	// https://de.wikipedia.org/wiki/Astronomische_Koordinatensysteme
	Equatorial_resting er;
	er.t = atan2(sin(hz.a)*cos(hz.h),sin(latitude)*cos(hz.a)*cos(hz.h)+cos(latitude)*sin(hz.h));
	er.d = asin(clamp( sin(latitude)*sin(hz.h) - cos(latitude)*cos(hz.h)*cos(hz.a) ,real(-1.0),real(1.0)));
	return er;
}


Horizontal horizontal(const Equatorial_resting& er,real latitude)
{
	using namespace std;
	// https://de.wikipedia.org/wiki/Astronomische_Koordinatensysteme
	Horizontal hz;
	hz.a = atan2(sin(er.t)*cos(er.d),sin(latitude)*cos(er.t)*cos(er.d)-cos(latitude)*sin(er.d));
	hz.h = asin(clamp( sin(latitude)*sin(er.d) + cos(latitude)*cos(er.d)*cos(er.t) ,real(-1.0),real(1.0)));
	return hz;
}




Ecliptic ecliptic(const Equatorial& eq, const real e)
{
	using namespace std;
	// https://de.wikipedia.org/wiki/Astronomische_Koordinatensysteme
	// https://en.wikipedia.org/wiki/Celestial_coordinate_system

	Ecliptic ec;
	ec.l = atan2(sin(eq.ra)*cos(e)*cos(eq.de) + sin(eq.de)*sin(e) ,cos(eq.ra)*cos(eq.de));
	ec.b = asin(clamp(cos(e)*sin(eq.de) - sin(e)*cos(eq.de)*sin(eq.ra),real(-1.0),real(1.0)));
	return ec;
}



Equatorial equatorial(const Ecliptic& ec,const real e)
{
	using namespace std;
	// https://de.wikipedia.org/wiki/Astronomische_Koordinatensysteme

	Equatorial eq;
	eq.ra = atan2(cos(e)*sin(ec.l)*cos(ec.b) - sin(e)*sin(ec.b), cos(ec.l)*cos(ec.b));
	eq.de = asin(clamp( cos(e)*sin(ec.b) + sin(e)*cos(ec.b)*sin(ec.l) ,real(-1.0),real(1.0)));
	return eq;
}



namespace {
	// https://en.wikipedia.org/wiki/Celestial_coordinate_system#Equatorial_%E2%86%94_galactic
	const real ra_p = rad(192.85948);
	const real de_p = rad( 27.12825);
	const real  l_p = rad(122.93192);
}



Equatorial equatorial(const Galactic& gc)
{
	using namespace std;
	
	Equatorial eq;
	eq.ra = ra_p + atan2(cos(gc.b)*sin(l_p-gc.l), cos(de_p)*sin(gc.b)-sin(de_p)*cos(gc.b)*cos(l_p-gc.l));
	eq.de = asin(clamp(sin(de_p)*sin(gc.b) + cos(de_p)*cos(gc.b)*cos(l_p-gc.l),real(-1.),real(1.)));
	return eq;
}



Galactic galactic(const Equatorial& eq)
{
	using namespace std;
	
	Galactic gc;
	gc.l = l_p - atan2(cos(eq.de)*sin(eq.ra-ra_p), cos(de_p)*sin(eq.de)-sin(de_p)*cos(eq.de)*cos(eq.ra-ra_p));
	gc.b = asin(clamp(sin(de_p)*sin(eq.de)+cos(de_p)*cos(eq.de)*cos(eq.ra-ra_p),real(-1.),real(1.)));
	return gc;
}



double sidereal_time(const Calendar& cal,double longitude)
{	
	// https://de.wikipedia.org/wiki/Sternzeit

	auto d0 = cal;
	auto UT = std::modf(cal.day,&d0.day);

	auto JD = julian_date(d0);
	auto T = julian_century(JD,Epoch_J2000);
	auto GMST = 100.46061837 + (36000.770053608 + (0.000387933 - 1.0/38710000.0*T)*T)*T; // (deg)
	auto t = GMST + (UT*360)*1.00273790935;

	return period_2pi(rad(t) + longitude);
}




Ecliptic sun_coordinates(double JD)
{
	// https://de.wikipedia.org/wiki/Sonnenstand

	auto n = JD - Epoch_J2000;
	auto L = 280.460 + 0.9856474*n; // (deg)
	auto g = 357.528 + 0.9856003*n; // (deg)

	auto l = L + 1.915*std::sin(rad(g)) + 0.01997*std::sin(rad(2*g)); // (deg)
	return {rad(period_360(l)), 0.0};
}




double obliquity_of_the_ecliptic(double JD)
{	
	// https://de.wikipedia.org/wiki/Ekliptik#Aktueller_Stand_der_Theorie
	
	double T = julian_century(JD,Epoch_J2000);
	double eps = 23.4392911111-(0.0130041667-(0.000000164+0.0000005036*T)*T)*T; // (deg)
	return rad(eps);
}



Equatorial earth_apex(double JD)
{
	auto s = sun_coordinates(JD);
	s.l -= rad(90);
	return equatorial(s,obliquity_of_the_ecliptic(JD));
}



}//aether
