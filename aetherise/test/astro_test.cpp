
#include "test.h"

#include "../astro.h"
#include "../DataSheet.h"

using namespace aether;


const double Obliquity_of_the_ecliptic = rad(23.44);

TEST(horizontal,equatorial_resting_to_horizontal)
{
	{
		Equatorial_resting ae;
		ae.t = rad(0);
		ae.d = rad(0);
		real latitude = rad(0);// equator
		Horizontal hz = horizontal(ae,latitude);
		//ASSERT_TRUE(test::approximate(hz.a,rad(0)));  // actually indeterminate
		ASSERT_TRUE(test::approximate(hz.h,rad(90)));
	}

	{
		Equatorial_resting ae;
		ae.t = rad(90);
		ae.d = rad(0);
		real latitude = rad(0);// equator
		Horizontal hz = horizontal(ae,latitude);
		ASSERT_TRUE(test::approximate(hz.a,rad(90)));
		ASSERT_TRUE(test::approximate(hz.h,rad(0)));
	}

	{
		Equatorial_resting ae;
		ae.t = rad(-90);
		ae.d = rad(0);
		real latitude = rad(0);// equator
		Horizontal hz = horizontal(ae,latitude);
		ASSERT_TRUE(test::approximate(hz.a,rad(-90)));
		ASSERT_TRUE(test::approximate(hz.h,rad(0)));
	}

	{
		Equatorial_resting ae;
		ae.t = rad(0);
		ae.d = rad(90);
		real latitude = rad(0);// equator
		Horizontal hz = horizontal(ae,latitude);
		ASSERT_TRUE(approximate_delta_periodic(hz.a,rad(180),AETHER_2PI,0.000001));
		ASSERT_TRUE(test::approximate(hz.h,rad(0)));
	}

	{
		Equatorial_resting ae;
		ae.t = rad(180);
		ae.d = rad(0);
		real latitude = rad(0);// equator
		Horizontal hz = horizontal(ae,latitude);
		//ASSERT_TRUE(test::approximate(hz.a,rad(90)));  // actually indeterminate
		ASSERT_TRUE(test::approximate(hz.h,rad(-90)));
	}

	{
		Equatorial_resting ae;
		ae.t = rad(-180);
		ae.d = rad(0);
		real latitude = rad(0);// equator
		Horizontal hz = horizontal(ae,latitude);
		//ASSERT_TRUE(test::approximate(hz.a,rad(-90)));  // actually indeterminate
		ASSERT_TRUE(test::approximate(hz.h,rad(-90)));
	}

	{
		Equatorial_resting ae;
		ae.t = rad(270);
		ae.d = rad(0);
		real latitude = rad(0);// equator
		Horizontal hz = horizontal(ae,latitude);
		ASSERT_TRUE(test::approximate(hz.a,rad(-90)));  // oder 270 normiert
		ASSERT_TRUE(test::approximate(hz.h,rad(0)));
	}

	{
		Equatorial_resting ae;
		ae.t = rad(0);
		ae.d = rad(0);
		real latitude = rad(45);// between north pole and equator
		Horizontal hz = horizontal(ae,latitude);
		ASSERT_TRUE(test::approximate(hz.a,rad(0)));
		ASSERT_TRUE(test::approximate(hz.h,rad(45)));
	}

	{
		Equatorial_resting ae;
		ae.t = rad(180);
		ae.d = rad(45);
		real latitude = rad(45);// between north pole and equator
		Horizontal hz = horizontal(ae,latitude);
		ASSERT_TRUE(approximate_delta_periodic(hz.a,rad(180),AETHER_2PI,0.000001));
		ASSERT_TRUE(test::approximate(hz.h,rad(0)));
	}

	{
		Equatorial_resting ae;
		ae.t = rad(0);
		ae.d = rad(-45);
		real latitude = rad(45);// between north pole and equator
		Horizontal hz = horizontal(ae,latitude);
		ASSERT_TRUE(test::approximate(hz.a,rad(0)));
		ASSERT_TRUE(test::approximate(hz.h,rad(0)));
	}

	{
		Equatorial_resting ae;
		ae.t = rad(0);
		ae.d = rad(0);
		real latitude = rad(90);// north pole
		Horizontal hz = horizontal(ae,latitude);
		ASSERT_TRUE(test::approximate(hz.a,rad(0)));
		ASSERT_TRUE(test::approximate(hz.h,rad(0)));
	}

	{
		Equatorial_resting ae;
		ae.t = rad(0);
		ae.d = rad(33);
		real latitude = rad(90);// north pole
		Horizontal hz = horizontal(ae,latitude);
		ASSERT_TRUE(test::approximate(hz.a,rad(0)));
		ASSERT_TRUE(test::approximate(hz.h,rad(33)));
	}


	{
		Equatorial_resting ae;
		ae.t = rad(0);
		ae.d = rad(-66);
		real latitude = rad(90);// north pole
		Horizontal hz = horizontal(ae,latitude);
		ASSERT_TRUE(test::approximate(hz.a,rad(0)));
		ASSERT_TRUE(test::approximate(hz.h,rad(-66)));
	}

	{
		Equatorial_resting ae;
		ae.t = rad(45);
		ae.d = rad(0);
		real latitude = rad(90);// north pole
		Horizontal hz = horizontal(ae,latitude);
		ASSERT_TRUE(test::approximate(hz.a,rad(45)));
		ASSERT_TRUE(test::approximate(hz.h,rad(0)));
	}

	{
		Equatorial_resting ae;
		ae.t = rad(90+45);
		ae.d = rad(0);
		real latitude = rad(90);// north pole
		Horizontal hz = horizontal(ae,latitude);
		ASSERT_TRUE(test::approximate(hz.a,rad(90+45)));
		ASSERT_TRUE(test::approximate(hz.h,rad(0)));
	}

	{
		Equatorial_resting ae;
		ae.t = rad(180+45);
		ae.d = rad(0);
		real latitude = rad(90);// north pole
		Horizontal hz = horizontal(ae,latitude);
		ASSERT_TRUE(test::approximate(hz.a,rad(-(180-45))));
		ASSERT_TRUE(test::approximate(hz.h,rad(0)));
	}

	{
		Equatorial_resting ae;
		ae.t = rad(-45);
		ae.d = rad(0);
		real latitude = rad(90);// north pole
		Horizontal hz = horizontal(ae,latitude);
		ASSERT_TRUE(test::approximate(hz.a,rad(-45)));
		ASSERT_TRUE(test::approximate(hz.h,rad(0)));
	}
}



TEST(rotating,test)
{
	{
		Equatorial_resting er = {1,2};
		Equatorial eq = rotating(er,3);
		ASSERT_APPROX(eq.ra,2,0.000001);
		ASSERT_APPROX(eq.de,2,0.000001);
	}
	
	{
		Equatorial_resting er = {1,2};
		Equatorial eq = rotating(er,-3);
		ASSERT_APPROX(eq.ra,-4,0.000001);
		ASSERT_APPROX(eq.de,2,0.000001);
	}
}



TEST(resting,test)
{
	{
		Equatorial eq {1,2};
		Equatorial_resting er = resting(eq,3);
		ASSERT_APPROX(er.t,2,0.00001);
		ASSERT_APPROX(er.d,2,0.00001);
	}
	
	
	{
		Equatorial eq {1,2};
		Equatorial_resting er = resting(eq,-3);
		ASSERT_APPROX(er.t,-4,0.00001);
		ASSERT_APPROX(er.d,2,0.00001);
	}
}



TEST(horizontal_equatorial,back_and_forth)
{
	{
		Horizontal hz;
		hz.a = rad(30);
		hz.h = rad(40);
		Equatorial_resting eqr = equatorial(hz,rad(50));
		Equatorial eq = rotating(eqr,h_to_rad(7));

		Horizontal hz2 = horizontal(resting(eq,h_to_rad(7)),rad(50));
		ASSERT_APPROX(hz2.a,hz.a,0.0000001);
		ASSERT_APPROX(hz2.h,hz.h,0.0000001);
	}

	{
		Equatorial eq;
		eq.ra = rad(30);
		eq.de = rad(40);
		Equatorial_resting eqr = resting(eq,h_to_rad(7));
		Horizontal hz = horizontal(eqr,rad(50));

		Equatorial eq2 = rotating(equatorial(hz,rad(50)),h_to_rad(7));
		ASSERT_APPROX(eq2.ra,eq.ra,0.0000001);
		ASSERT_APPROX(eq2.de,eq.de,0.0000001);
	}
}






TEST(equatorial,ecliptic)
{
	const double oote = rad(23.44);
	{		
		Ecliptic ek {rad(0),rad(0)};
		Equatorial ae = equatorial(ek,oote);
		ASSERT_APPROX(ae.ra,rad(0),0.000001);
		ASSERT_APPROX(ae.de,rad(0),0.000001);
	}
	
	{		
		Ecliptic ek {rad(180),rad(0)};
		Equatorial ae = equatorial(ek,oote);
		ASSERT_APPROX(period_2pi(ae.ra),rad(180),0.000001);
		ASSERT_APPROX(ae.de,rad(0),0.000001);
	}
	
	{		
		Ecliptic ek {rad(90),rad(0)};
		Equatorial ae = equatorial(ek,oote);
		ASSERT_APPROX(ae.ra,rad(90),0.000001);
		ASSERT_APPROX(ae.de,oote,0.000001);
	}
	
	{		
		Ecliptic ek {rad(-90),rad(0)};
		Equatorial ae = equatorial(ek,oote);
		ASSERT_APPROX(period_2pi(ae.ra),rad(360-90),0.000001);
		ASSERT_APPROX(ae.de,-oote,0.000001);
	}
}



TEST(ecliptic,test)
{
	{
		Equatorial ae {h_to_rad(0),rad(0)};
		Ecliptic ek = ecliptic(ae,Obliquity_of_the_ecliptic);
		ASSERT_APPROX(ek.l,rad(0),0.000001);
		ASSERT_APPROX(ek.b,rad(0),0.000001);
	}

	{
		Equatorial ae {h_to_rad(12),rad(0)};
		Ecliptic ek = ecliptic(ae,Obliquity_of_the_ecliptic);
		ASSERT_APPROX(ek.l,rad(180),0.000001);
		ASSERT_APPROX(ek.b,rad(0)  ,0.000001);
	}

	{
		Equatorial ae {h_to_rad(6),rad(0)};
		Ecliptic ek = ecliptic(ae,Obliquity_of_the_ecliptic);
		ASSERT_APPROX(ek.l,rad(90),0.000001);
		ASSERT_APPROX(ek.b,-Obliquity_of_the_ecliptic,0.000001);
	}
	
	{
		Equatorial ae {h_to_rad(-6),rad(0)};
		Ecliptic ek = ecliptic(ae,Obliquity_of_the_ecliptic);
		ASSERT_APPROX(period_2pi(ek.l),rad(360-90),0.000001);
		ASSERT_APPROX(ek.b,Obliquity_of_the_ecliptic,0.000001);
	}
}



TEST(ecliptic,equatorial_ecliptic_transformation)
{
	{
		for (int ra=0;ra<=360;ra+=15) {
			for (int de=-75;de<=75;de+=15) {
				Equatorial ae {rad(ra),rad(de)};
				Ecliptic ek = ecliptic(ae,Obliquity_of_the_ecliptic);
				auto ae2 = equatorial(ek,Obliquity_of_the_ecliptic);
				ASSERT_TRUE(approximate_delta_periodic(ae2.ra,ae.ra,0.000000001,AETHER_2PI));
				ASSERT_TRUE(approximate(ae2.de,ae.de,0.000000001));
			}
		}
	}
}



TEST(sidereal_time,test)
{
	{
		Calendar cal {1925,9,23 + (3. + 3./60)/24};
		auto theta = sidereal_time(cal,0);
		auto h = theta/AETHER_PI*12;
		// at 3:09 its ~ 3h 22m
		ASSERT_TRUE(approximate_delta(h,3,0.25));
	}
	
	{
		Calendar cal {1925,9,23 + (3. +7.86 + 3./60)/24};
		auto theta = sidereal_time(cal,rad(-118));
		auto h = theta/AETHER_PI*12;
		// at local 3:09 on lon -118 its ~ 3h 22m
		ASSERT_TRUE(approximate_delta(h,3,0.25));
	}
	
	{
		// more tests via -validate with aetherise
	}
}



TEST(conversion,equatorial_polar)
{
	{
		double r1 = 1;
		double ra1 = 0;
		double de1 = 0;

		double r2 = 1;
		double ra2 = 0;
		double de2 = 0;

		Polar3 p1;
		p1.r = r1;
		p1.th = rad(90)-de1;
		p1.ph = ra1;

		Polar3 p2;
		p2.r = r2;
		p2.th = rad(90)-de2;
		p2.ph = ra2;

		auto v1 = vector3(p1);
		auto v2 = vector3(p2);
		auto v = v1+v2;
		auto p = polar3(v);

		double r3 = p.r;
		double ra3 = p.ph;
		double de3 = rad(90)-p.th;

		ASSERT_APPROX(r3,2,1e-6);
		ASSERT_APPROX(ra3,0,1e-6);
		ASSERT_APPROX(de3,0,1e-6);
	}

	{
		double r1 = 1;
		double ra1 = 0;
		double de1 = 0;

		double r2 = 1;
		double ra2 = h_to_rad(6);
		double de2 = 0;

		Polar3 p1;
		p1.r = r1;
		p1.th = rad(90)-de1;
		p1.ph = ra1;

		Polar3 p2;
		p2.r = r2;
		p2.th = rad(90)-de2;
		p2.ph = ra2;

		auto v1 = vector3(p1);
		auto v2 = vector3(p2);
		auto v = v1+v2;
		auto p = polar3(v);

		double r3 = p.r;
		double ra3 = p.ph;
		double de3 = rad(90)-p.th;

		ASSERT_APPROX(r3,std::sqrt(2),1e-6);
		ASSERT_APPROX(ra3,h_to_rad(3),1e-6);
		ASSERT_APPROX(de3,0,1e-6);
	}


	{
		double r1 = 1;
		double ra1 = h_to_rad(-6);
		double de1 = rad(45);

		double r2 = 1;
		double ra2 = h_to_rad(18);
		double de2 = rad(-45);

		Polar3 p1;
		p1.r = r1;
		p1.th = rad(90)-de1;
		p1.ph = ra1;

		Polar3 p2;
		p2.r = r2;
		p2.th = rad(90)-de2;
		p2.ph = ra2;

		auto v1 = vector3(p1);
		auto v2 = vector3(p2);
		auto v = v1+v2;

		ASSERT_APPROX(length(v1),1,1e-6);
		ASSERT_APPROX(length(v2),1,1e-6);
		ASSERT_APPROX(length(v),std::sqrt(2),1e-6);

		auto p = polar3(v);
		double r3 = p.r;
		double ra3 = p.ph;
		double de3 = rad(90)-p.th;

		ASSERT_APPROX(r3,std::sqrt(2),1e-6);
		ASSERT_APPROX(period_2pi(ra3),h_to_rad(18),1e-6);
		ASSERT_APPROX(de3,0,1e-6);
	}
}



TEST(obliquity_of_the_ecliptic,test)
{
	{
		auto JD = julian_date({2021,5,25});
		auto eps = obliquity_of_the_ecliptic(JD);
		ASSERT_APPROX(eps,rad(23.43651),0.000001); // wikipedia 25.5.2021
	}
}



TEST(sun_coordinates,test)
{
		
	// values from calsky.com
	{
		// timezone & daylight saving time needs to be recognised?
		Calendar cal {2020,5,25+(16-2)/24.0}; // 25.5.2020 16:00 in calsky.com
		auto JD = julian_date(cal);
		auto ek = sun_coordinates(JD);

		real expected_l = rad(64 + 49/60.0 + 33.4/3600.0);
		ASSERT_TRUE(approximate(ek.l,expected_l,0.001)); // accuracy?
	}


	{
		Calendar cal {1925,7,30+((5-2)+10/60.0)/24.0}; // 30.7.1925 5:10 in calsky.com
		auto JD = julian_date(cal);
		//JD += 7.86/24;
		auto ek = sun_coordinates(JD);

		real expected_l = rad(126 + 28/60.0 + 37/3600.0);
		ASSERT_TRUE(approximate(ek.l,expected_l,0.001)); // accuracy?
	}
	
	{
		Calendar cal {1925,7,30+((5)+10/60.0)/24.0}; // 30.7.1925 5:10 in UT
		auto JD = julian_date(cal);
		JD += 7.86/24;
		auto ek = sun_coordinates(JD);
		auto ae = equatorial(ek,obliquity_of_the_ecliptic(JD));
		auto theta = h_to_rad(time_to_h({1,37}));
		auto hz = horizontal(resting(ae,theta),rad(34.225));
		
		std::cout << "Mt. Wilson, sunrise\n";
		std::cout << "Azimuth: " << deg(hz.a) << " (south)\n";
		std::cout << "Height: " << deg(hz.h) << "\n";
	}
	
	
	{
		// https://aa.quae.nl/en/reken/zonpositie.html
		
		Calendar cal {2004,4,1+(12/24.0)}; // 1.4.2004 12:00 UTC
		auto JD = julian_date(cal);
		auto ek = sun_coordinates(JD);
		auto eq = equatorial(ek,Obliquity_of_the_ecliptic);		
		auto theta = sidereal_time(cal,rad(5));
		auto hz = horizontal(resting(eq,theta),rad(52));

		real expected_a = rad(5);
		real expected_h = rad(43);
		ASSERT_TRUE(approximate(hz.a,expected_a,0.2)); // not very accurate
		ASSERT_TRUE(approximate(hz.h,expected_h,0.2)); 
	}
	

	
	{
		// https://www.sunearthtools.com/dp/tools/pos_sun.php?lang=de
		
		double lat = 52.0398774;
		double lon = 6.9854455;
		Calendar cal {2004,4,1+(12/24.0)}; // 
		auto JD = julian_date(cal);
		auto ek = sun_coordinates(JD);
		auto eq = equatorial(ek,Obliquity_of_the_ecliptic);		
		auto theta = sidereal_time(cal,rad(lon));
		auto hz = horizontal(resting(eq,theta),rad(lat));

		real expected_a = rad(188.17);
		real expected_h = rad(42.48);
		auto na = north_azimuth(hz.a);
		ASSERT_TRUE(approximate(na,expected_a,0.01)); // better accuracy
		ASSERT_TRUE(approximate(hz.h,expected_h,0.01)); 
	}
	
}



TEST(earth_apex,test)
{
	{
		// winter solstice
		Calendar cal {2019,12,21+(12)/24.0}; // 
		auto JD = julian_date(cal);
		Equatorial eq = earth_apex(JD);
		auto theta = sidereal_time(cal,rad(0));
		auto hz = horizontal(resting(eq,theta),rad(66.6));
				
		ASSERT_APPROX(period_2pi(hz.a),rad(90),0.01);
		ASSERT_APPROX(hz.h,rad(0),0.01);		
	}
	
	{
		// Mt. Wilson, sunrise
		Calendar cal {1925,7,30+((5+7.86)+10/60.0)/24.0}; // 30.7.1925 5:10 local time in UT
		auto JD = julian_date(cal);
		Equatorial eq = earth_apex(JD);
		auto theta = h_to_rad(time_to_h({1,37}));
		auto hz = horizontal(resting(eq,theta),rad(34.225));
				
		ASSERT_TRUE(approximate_delta(period_2pi(hz.a),rad(360-118+90),rad(5)));
		ASSERT_TRUE(hz.h > rad(45));		
	}
}



int main(int ,char** )
{
	RUN_ALL_TESTS();

	return EXIT_SUCCESS;
}
