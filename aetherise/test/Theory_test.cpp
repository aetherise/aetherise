
#include "test.h"

#include "../Theory.h"
#include "../astro.h"

using namespace aether;

#define MICRO 0.000001
#define NANO 0.000000001

TEST(Classic,test)
{

	{
		Classic theorie;
		const auto c = Lightspeed;

		real V = 300000;
		real L = 1;
		Vector3 v {V,0,0};

		Lightpath lp;		
		lp.n = 1;
		lp.p = Vector3 {L,0,0}; // parallel to movement

		real t = theorie.phase_delay(v,lp);
		real expected_t = L/(c-V);
		ASSERT_TRUE(approximate(t,expected_t,NANO));

	}

	{
		Classic theorie;
		const auto c = Lightspeed;

		real V = 300000;
		real L = 1;
		Vector3 v {V,0,0};

		Lightpath lp;		
		lp.n = 1;
		lp.p = Vector3 {-L,0,0}; // parallel to movement

		real t = theorie.phase_delay(v,lp);
		real expected_t = L/(c+V);
		ASSERT_TRUE(approximate(t,expected_t,NANO));
	}


	{
		Classic theorie;
		const auto c = Lightspeed;

		real V = 300000;
		real L = 1;
		Vector3 v {V,0,0};

		Lightpath lp;		
		lp.n = 1;
		lp.p = Vector3 {0,L,0}; // perpendicular to movement

		real t = theorie.phase_delay(v,lp);
		real expected_t = L/c * 1/sqrt_1_minus_sqr(V/c);
		ASSERT_TRUE(approximate(t,expected_t,NANO));
	}

	{
		Classic theorie;
		const auto c = Lightspeed;

		real V = 300000;
		real L = 1;
		Vector3 v {V,0,0};

		Lightpath lp;		
		lp.n = 1;
		lp.p = Vector3 {0,-L,0}; // perpendicular to movement

		real t = theorie.phase_delay(v,lp);
		real expected_t = L/c * 1/sqrt_1_minus_sqr(V/c);
		ASSERT_TRUE(approximate(t,expected_t,NANO));
	}
}






TEST(Aether,test)
{
	{
		Aether aether;

		real V = 300000;
		real L = 1;
		double n = 1; // Vacuum
		Vector3 v {V,0,0};

		std::vector<Lightpath> lp1 {{{0,L,0},n}, {{0,-L,0},n}};
		std::vector<Lightpath> lp2 {{{L,0,0},n}, {{-L,0,0},n}};

		real t1 = phase_delay(aether,v,lp1);
		real t2 = phase_delay(aether,v,lp2);

		ASSERT_APPROX(t1,t2,NANO);
	}

	{
		Aether aether;

		real V = 300000;
		real L = 1;
		double n = 1; // Vacuum
		Vector3 v {0,V,0};

		std::vector<Lightpath> lp1 {{{0,L,0},n}, {{0,-L,0},n}};
		std::vector<Lightpath> lp2 {{{L,0,0},n}, {{-L,0,0},n}};

		real t1 = phase_delay(aether,v,lp1);
		real t2 = phase_delay(aether,v,lp2);

		ASSERT_APPROX(t1,t2,NANO);
	}

	{
		Aether aether;

		real V = 300000;
		real L = 1;
		double n = 1.5; // impossible gas
		Vector3 v {0,0,V}; // no component of v in plane of interferometer

		std::vector<Lightpath> lp1 {{{0,L,0},n}, {{0,-L,0},n}};
		std::vector<Lightpath> lp2 {{{L,0,0},n}, {{-L,0,0},n}};

		real t1 = phase_delay(aether,v,lp1);
		real t2 = phase_delay(aether,v,lp2);

		ASSERT_APPROX(t1,t2,NANO);
	}

	{
		Aether aether;

		real V = 300000;
		real L = 1;
		double n = 1.5; // impossible gas
		Vector3 v {V,0,0};

		std::vector<Lightpath> lp1 {{{0,L,0},n}, {{0,-L,0},n}};
		std::vector<Lightpath> lp2 {{{L,0,0},n}, {{-L,0,0},n}};

		real t1 = phase_delay(aether,v,lp1);
		real t2 = phase_delay(aether,v,lp2);

		ASSERT_TRUE(!approximate(t1,t2,NANO));
	}
}




TEST(phase_delay,test)
{
	{
		Classic classic;
		Vector3 v {300000,0,0};
		std::vector<Lightpath> lp;
		auto t = phase_delay(classic,v,lp);
		ASSERT_TRUE(t==0);
	}

	{
		Classic classic;
		Vector3 v {300000,0,0};
		real L = 1;
		std::vector<Lightpath> lp {{{0,L,0},1}, {{0,-L,0},1}};
		auto t = phase_delay(classic,v,lp);
		auto s = classic.phase_delay(v,lp.at(0)) + classic.phase_delay(v,lp.at(1));
		ASSERT_TRUE(t==s);
	}
}




TEST(fringe_displacements, sign_correct)
{
	// According to Miller, the note 'sign correct' means
	// that lengthening the telescope arm gives increasing positive readings.
	// So we test our calculations here and set the correct sign
	// in fringe_displacements().
	{
		const real L = 32;
		const double n = 1; // vacuum
		const Vector3 v = {0,370000,0};
		Aether theory;

		real dt1,dt2;
		{
			Interferometer interferometer;
			interferometer.light_path_1 = { {{0,L,0},n}, {{0,-L,0},n} };
			interferometer.light_path_2 = { {{L,0,0},n}, {{-L,0,0},n} };
			interferometer.wave_length = 570e-9;

			real t1 = phase_delay(theory,v,interferometer.light_path_1);
			real t2 = phase_delay(theory,v,interferometer.light_path_2);
			dt1 = -(t2-t1);
		}

		{
			real L2 = L*1.000001;
			Interferometer interferometer;
			interferometer.light_path_1 = { {{0,L2,0},n}, {{0,-L2,0},n} };
			interferometer.light_path_2 = { {{L,0,0},n}, {{-L,0,0},n} };
			interferometer.wave_length = 570e-9;

			real t1 = phase_delay(theory,v,interferometer.light_path_1);
			real t2 = phase_delay(theory,v,interferometer.light_path_2);
			dt2 = -(t2-t1);
		}

		ASSERT_TRUE(dt2-dt1 > 0);
	}
}




TEST(fringe_displacements,relativity)
{	
	const double lat = MtWilson_Latitude;
	const double L = Millers_Interferometer_Arm_Length;
	
	{
		Relativity theory;
		std::array<double,17> displs = fringe_displacements(theory,CMB_dipole,lat,1.00023,L,0,false);
		for (auto displ : displs) {
			ASSERT_APPROX(displ,0,NANO);
		}
	}
}



TEST(fringe_displacements,ether)
{
	const double lat = MtWilson_Latitude;
	const double L = Millers_Interferometer_Arm_Length;
	
	{
		// On summer solstice, a movement in direction of the sun at noon
		// should result in an almost zero amplitude of the signal.
		Aether theory;
		
		Calendar cal {1925,6,21 + (12.)/24.};				
		auto JD = julian_date(cal);
		auto ec = sun_coordinates(JD);
		auto eq = equatorial(ec,rad(23.44));
		TheoryParameters tpsun {370000,double(eq.ra),double(eq.de)};
		auto theta = sidereal_time(cal,0);
		std::array<double,17> displs = fringe_displacements(theory,tpsun,lat,1.00023,L,theta,false);
		auto amp = max_abs_value(displs);
		
		ASSERT_TRUE(amp < 0.005);
	}	
	
	{
		// On summer solstice, a movement perpendicular to the direction of the sun at noon
		// should result in an almost maximal amplitude of the signal.
		Aether theory;
		
		Calendar cal {1925,6,21 + (12.)/24.};				
		auto JD = julian_date(cal);
		auto ec = sun_coordinates(JD);
		ec.b = rad(90); // make perpendicular
		auto eq = equatorial(ec,rad(23.44));
		TheoryParameters tpsun {370000,double(eq.ra),double(eq.de)};
		auto theta = sidereal_time(cal,0);
		std::array<double,17> displs = fringe_displacements(theory,tpsun,lat,1.00023,L,theta,false);
		auto amp = max_abs_value(displs);
		
		ASSERT_TRUE(amp > 0.015);
		// Because of Miller's definition of 'sign correct', the first values
		// are expected to be negative.
		ASSERT_TRUE(displs.at(0)<0 && displs.at(1)<0);
	}	
}



int main(int ,char** )
{
	RUN_ALL_TESTS();

	return EXIT_SUCCESS;
}
