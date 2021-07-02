#include "test.h"

#include "../aetherise.h"
#include "../mathematics.h"
#include "../mini.h"
#include "../astro.h"
#include "../Theory.h"

#include <random>
#include <iostream>
#include <iomanip>

using namespace aether;


TEST(subtract_data,test)
{
	{
		Options options;
		options.subtract_data = true;
		options.data1 = {2,2,2,2, 2,2,2,2, 2,2,2,2, 2,2,2,2, 2};
		options.data2 = {.1,.1,.1,.1, .1,.1,.1,.1, .1,.1,.1,.1, .1,.1,.1,.1, .1};
		options.data3 = {1,1,1,1, 1,1,1,1, 1,1,1,1, 1,1,1,1, 1};
					
		DisplacementData displs;
		displs.data = {1,1,1,1, 1,1,1,1, 1,1,1,1, 1,1,1,1, 1};
		displs.uncertainties = {.1,.1,.1,.1, .1,.1,.1,.1, .1,.1,.1,.1, .1,.1,.1,.1, .1};			   
		displs.theory = {1,1,1,1, 1,1,1,1, 1,1,1,1, 1,1,1,1, 1};
		
		subtract_data(displs,options);
		
		for (auto& d : displs.data) {
			ASSERT_TRUE(d==-1);
		}
		for (auto& u : displs.uncertainties) {
			ASSERT_APPROX(u,std::sqrt(2*0.1*0.1),0.000001);
		}
		for (auto& t : displs.theory) {
			ASSERT_TRUE(t==0);
		}
	}
}



TEST(subtract_data,exit)
{
	{
		Options options;
		options.subtract_data = true;
		options.data1 = {2,2,2,2, 2,2,2,2, 2,2,2,2, 2,2,2,2};// missing
		options.data2 = {.1,.1,.1,.1, .1,.1,.1,.1, .1,.1,.1,.1, .1,.1,.1,.1 ,.1}; 
		options.data3 = {1,1,1,1, 1,1,1,1, 1,1,1,1, 1,1,1,1, 1};
					
		DisplacementData displs;
				
		try {
			subtract_data(displs,options);
			FAIL();
		}
		catch(ExitException&) {			
		}		
	}
	
	{
		Options options;
		options.subtract_data = true;
		options.data1 = {2,2,2,2, 2,2,2,2, 2,2,2,2, 2,2,2,2, 2};
		options.data2 = {.1,.1,.1,.1, .1,.1,.1,.1, .1,.1,.1,.1, .1,.1,.1,.1}; // missing
		options.data3 = {1,1,1,1, 1,1,1,1, 1,1,1,1, 1,1,1,1, 1};
					
		DisplacementData displs;
				
		try {
			subtract_data(displs,options);
			FAIL();
		}
		catch(ExitException&) {			
		}		
	}
	
	{
		Options options;
		options.subtract_data = true;
		options.data1 = {2,2,2,2, 2,2,2,2, 2,2,2,2, 2,2,2,2, 2};
		options.data2 = {.1,.1,.1,.1, .1,.1,.1,.1, .1,.1,.1,.1, .1,.1,.1,.1, .1}; 
		options.data3 = {1,1,1,1, 1,1,1,1, 1,1,1,1, 1,1,1,1}; // missing
					
		DisplacementData displs;
				
		try {
			subtract_data(displs,options);
			FAIL();
		}
		catch(ExitException&) {			
		}		
	}
	
	{
		Options options;
		options.subtract_data = true;
		options.model = true;
		options.data1 = {2,2,2,2, 2,2,2,2, 2,2,2,2, 2,2,2,2, 2};
		options.data2 = {.1,.1,.1,.1, .1,.1,.1,.1, .1,.1,.1,.1, .1,.1,.1,.1, .1}; 
		options.data3 = {1,1,1,1, 1,1,1,1, 1,1,1,1, 1,1,1,1}; 
		// missing model
					
		DisplacementData displs;
				
		try {
			subtract_data(displs,options);
			FAIL();
		}
		catch(ExitException&) {			
		}		
	}
}



TEST(add_earth_orbit,test)
{
	const Options options;
	
	{	
		// Try to cancel out earth velocity
		// by assuming movement in the ether in the opposite direction.
		DataSheet data_sheet;
		data_sheet.date = {1925,6,21};
		data_sheet.mean_observation_time = {12+8,0};
		Calendar cal = calendar_date(data_sheet);
		Horizontal hz {rad(-90),0};
		auto theta = sidereal_time(cal,rad(-118)); // Mt. Wilson
		Equatorial eq = rotating(equatorial(hz,rad(34)),theta);
		
		TheoryParameters params {29900,double(eq.ra),double(eq.de)}; // approx opposite apex
		
		auto weo = add_earth_orbit(params,data_sheet,options);
		ASSERT_TRUE(weo.v>0 && weo.v < 1000);
	}
	
	{	
		// Try to double result velocity
		// by assuming movement in the ether in the same direction.
		DataSheet data_sheet;
		data_sheet.date = {1925,6,21};
		data_sheet.mean_observation_time = {12+8,0};
		Calendar cal = calendar_date(data_sheet);
		Horizontal hz {rad(90),0};
		auto theta = sidereal_time(cal,rad(-118)); // Mt. Wilson
		Equatorial eq = rotating(equatorial(hz,rad(34)),theta);
		
		TheoryParameters params {29900,double(eq.ra),double(eq.de)}; // approx same apex
		
		auto weo = add_earth_orbit(params,data_sheet,options);
		ASSERT_APPROX(weo.v,60000,0.1);
	}
}



TEST(humidity,relative_validity)
{
	{
		DataSheet data_sheet;		
		auto h_clear = humidity(data_sheet);
		
		data_sheet.weather_obscuration = DataSheet::WeatherObscurations::Haze;
		auto h_haze = humidity(data_sheet);
		
		data_sheet.weather_obscuration = DataSheet::WeatherObscurations::Mist;
		auto h_mist = humidity(data_sheet);
		
		data_sheet.weather_obscuration = DataSheet::WeatherObscurations::Fog;
		auto h_fog = humidity(data_sheet);
		
		ASSERT_TRUE(h_clear <= h_haze);
		ASSERT_TRUE(h_haze < h_mist);
		ASSERT_TRUE(h_mist < h_fog);
		
		data_sheet.weather_rain = true;
		data_sheet.weather_obscuration = DataSheet::WeatherObscurations::Clear;
		auto h_clear_rain = humidity(data_sheet);
		ASSERT_TRUE(h_clear_rain > h_clear);
	}		
}



TEST(index_of_refraction,test)
{
	{
		DataSheet data_sheet;
		
		auto n = index_of_refraction(data_sheet);
		ASSERT_TRUE(n==1.00023);
	}
	
	{
		DataSheet data_sheet;
		data_sheet.thermometers_start = {15,15,15,15};
		auto n = index_of_refraction(data_sheet);
		ASSERT_TRUE(n!=1.00023);	
		ASSERT_APPROX(n,1.00023,0.00001);
		
		data_sheet.thermometers_start = {5,5,5,5};
		auto n2 = index_of_refraction(data_sheet);
		ASSERT_TRUE(n < n2);	
	}
	
	{
		DataSheet data_sheet;
		data_sheet.thermometers_end = {15,15,15,15};
		auto n = index_of_refraction(data_sheet);
		ASSERT_TRUE(n!=1.00023);	
		ASSERT_APPROX(n,1.00023,0.00001);
		
		data_sheet.thermometers_end = {5,5,5,5};
		auto n2 = index_of_refraction(data_sheet);
		ASSERT_TRUE(n < n2);	
	}
}



TEST(fringe_displacements,earth_enabled)
{
	{
		// At the end of august, the sun stands in the zodiak lion.
		// Thats approximate the apex of movement according to the CMB dipole.
		Aether theory;
		TheoryParameters params = CMB_dipole;
		DataSheet data_sheet;
		data_sheet.date = {1925,9-3,1};
		data_sheet.mean_observation_time = {12+6,0}; 
		Calendar cal = calendar_date(data_sheet);
		auto theta1 = sidereal_time(cal,0);
		data_sheet.sidereal_time = h_to_time(theta1);
		data_sheet.thermometers_start = {15,15,15,15};
		Options options;		
		auto displs1 = fringe_displacements(theory,params,data_sheet,options);
				
		data_sheet.date = {1925,9+3,1};
		data_sheet.mean_observation_time = {12-6,0}; 
		cal = calendar_date(data_sheet);
		auto theta2 = sidereal_time(cal,0);	
		
		ASSERT_APPROX(theta1,theta2,0.01);
		
		data_sheet.sidereal_time = h_to_time(theta2);
		auto displs2 = fringe_displacements(theory,params,data_sheet,options);
		
		auto amp1 = max_abs_value(displs1);
		auto amp2 = max_abs_value(displs2);
		// In summer the earth has an opposite apex, so the amplitude is lower.
		// In winter the earth has the same apex, so the amplitude is higher.		
		ASSERT_TRUE(amp1 < amp2);
	}
}



TEST(fringe_displacements,earth_disabled)
{
	{
		// At the end of august, the sun stands in the zodiak lion.
		// Thats approximate the apex of movement according to the CMB dipole.
		Aether theory;
		TheoryParameters params = CMB_dipole;
		DataSheet data_sheet;
		data_sheet.date = {1925,9-3,1};
		data_sheet.mean_observation_time = {12+6,0}; 
		Calendar cal = calendar_date(data_sheet);
		auto theta1 = sidereal_time(cal,0);
		data_sheet.sidereal_time = h_to_time(theta1);
		data_sheet.thermometers_start = {15,15,15,15};
		Options options;	
		options.enable_earth = false; // earth disabled
		auto displs1 = fringe_displacements(theory,params,data_sheet,options);
				
		data_sheet.date = {1925,9+3,1};
		data_sheet.mean_observation_time = {12-6,0}; 
		cal = calendar_date(data_sheet);
		auto theta2 = sidereal_time(cal,0);	
		
		ASSERT_APPROX(theta1,theta2,0.01);
		
		data_sheet.sidereal_time = h_to_time(theta2);
		auto displs2 = fringe_displacements(theory,params,data_sheet,options);
		
		auto amp1 = max_abs_value(displs1);
		auto amp2 = max_abs_value(displs2);
		// In summer the earth has an opposite apex, so the amplitude is lower.
		// In winter the earth has the same apex, so the amplitude is higher.		
		ASSERT_APPROX(amp1,amp2,0.01);
	}
}



TEST(rejection_quota,test)
{
	{
		auto ci1 = rejection_quota(25,100);
		auto ci2 = confidence_interval(25,100);
		ASSERT_APPROX(ci1.p,ci2.p,0.000001);
		ASSERT_APPROX(ci1.u,ci2.u,0.000001);
	}
}



TEST(estimate_amplitude,test)
{
	{
		std::array<double,17> data {0,0,0,1, 0,0,0,0, -3,0,0,0, 1,1,1,1, 0};
		auto amp = estimate_amplitude(data);
		ASSERT_APPROX(amp,2,0.1);
	}
}



TEST(sidereal_time_to_key,test)
{
	{
		auto key = sidereal_time_to_key(0);
		ASSERT_APPROX(key,0,0.000001);			
	}
	
	{
		auto key1 = sidereal_time_to_key(0);
		auto key2 = sidereal_time_to_key(0+Sidereal_Aggregation_Bin_Width);
		ASSERT_TRUE(key1+1 == key2);			
	}
	
	{
		auto key = sidereal_time_to_key(7);
		ASSERT_APPROX(key,7/Sidereal_Aggregation_Bin_Width,0.000001);			
	}
	
	{
		auto key1 = sidereal_time_to_key(7);
		auto key2 = sidereal_time_to_key(7+Sidereal_Aggregation_Bin_Width);
		ASSERT_TRUE(key1+1 == key2);			
	}
}



TEST(sidereal_time_to_key,key_to_mean_sidereal_time)
{
	{
		for (int i=0;i<24*5;i++) {
			double h = i/5.;
			ASSERT_TRUE(h>=0 && h<24);
			auto key = sidereal_time_to_key(h);
			ASSERT_TRUE(key>=0);
			auto h2 = key_to_mean_sidereal_time(key);
			ASSERT_TRUE(h>=h2-Sidereal_Aggregation_Bin_Width/2 && h<=h2+Sidereal_Aggregation_Bin_Width/2);
		}
	}
}



TEST(key_to_mean_sidereal_time,test)
{
	{
		auto h = key_to_mean_sidereal_time(0);
		ASSERT_APPROX(h,Sidereal_Aggregation_Bin_Width/2,0.00001);
	}
}





void fill_with_sine(double p, double a,double c,std::array<double,17>& data, std::array<double,17>& uncertainties)
{
	for (size_t i=0;i<17;i++) {
		data.at(i) = a*std::sin(i*AETHER_PI/4-p) + c;
		uncertainties.at(i) = 0.015;
	}
}



TEST(fit_sine,grad)
{
	Options options;
	options.single = true;
	options.minimizer = Options::Minimizer::Grad;
	
	{
		std::array<double,17> data;
		std::array<double,17> uncertainties;
		double p = rad(0);
		double a = 0;
		fill_with_sine(p,a,0,data,uncertainties);
		auto result = fit_sine(data,uncertainties,options);
		
		ASSERT_TRUE(result.valid);
		//ASSERT_APPROX(result.x[0],p,0.001);
		ASSERT_APPROX(result.x[1],a,0.001);
	}
	
	{
		std::array<double,17> data;
		std::array<double,17> uncertainties;
		double p = rad(0);
		double a = 0.02;
		fill_with_sine(p,a,0,data,uncertainties);
		auto result = fit_sine(data,uncertainties,options);
		
		ASSERT_TRUE(result.valid);
		ASSERT_APPROX(result.x[0],p,0.001);
		ASSERT_APPROX(result.x[1],a,0.001);
	}
	
	{
		std::array<double,17> data;
		std::array<double,17> uncertainties;
		double p = rad(45);
		double a = 0.02;
		fill_with_sine(p,a,0,data,uncertainties);
		auto result = fit_sine(data,uncertainties,options);
		
		ASSERT_TRUE(result.valid);
		ASSERT_APPROX(result.x[0],p,0.001);
		ASSERT_APPROX(result.x[1],a,0.001);
	}
	
	{
		std::array<double,17> data;
		std::array<double,17> uncertainties;
		double p = rad(45);
		double a = 0.04; // higher amp
		fill_with_sine(p,a,0,data,uncertainties);
		auto result = fit_sine(data,uncertainties,options);
		
		ASSERT_TRUE(result.valid);
		ASSERT_APPROX(result.x[0],p,0.001);
		ASSERT_APPROX(result.x[1],a,0.001);
	}
	
	{
		std::array<double,17> data;
		std::array<double,17> uncertainties;
		double p = rad(0);
		double a = 0.02;
		fill_with_sine(p,a,2.4,data,uncertainties);
		auto result = fit_sine(data,uncertainties,options);
		
		ASSERT_TRUE(result.valid);
		ASSERT_APPROX(result.x[0],p,0.001);
		ASSERT_APPROX(result.x[1],a,0.001);
	}
	
	{
		std::array<double,17> data;
		std::array<double,17> uncertainties;
		double p = rad(45);
		double a = 0.02;
		fill_with_sine(p,a,-2.4,data,uncertainties);
		auto result = fit_sine(data,uncertainties,options);
		
		ASSERT_TRUE(result.valid);
		ASSERT_APPROX(result.x[0],p,0.001);
		ASSERT_APPROX(result.x[1],a,0.001);
	}
	
	{
		std::array<double,17> data;
		std::array<double,17> uncertainties;
		double p = rad(-45);
		double a = 0.02;
		fill_with_sine(p,a,-2.4,data,uncertainties);
		auto result = fit_sine(data,uncertainties,options);
		
		ASSERT_TRUE(result.valid);
		ASSERT_APPROX(result.x[0],period_2pi(p),0.001);
		ASSERT_APPROX(result.x[1],a,0.001);
	}
	
	{
		std::array<double,17> data;
		std::array<double,17> uncertainties;
		double p = rad(180);
		double a = 0.02;
		fill_with_sine(p,a,-2.4,data,uncertainties);
		auto result = fit_sine(data,uncertainties,options);
		
		ASSERT_TRUE(result.valid);
		ASSERT_APPROX(result.x[0],p,0.001);
		ASSERT_APPROX(result.x[1],a,0.001);
	}
	
	{
		std::array<double,17> data;
		std::array<double,17> uncertainties;
		double p = rad(250);
		double a = 0.02;
		fill_with_sine(p,a,-2.4,data,uncertainties);
		auto result = fit_sine(data,uncertainties,options);
		
		ASSERT_TRUE(result.valid);
		ASSERT_APPROX(result.x[0],p,0.001);
		ASSERT_APPROX(result.x[1],a,0.001);
	}
}



#ifdef AETHER_MINUIT

TEST(fit_sine,Minuit2)
{
	Options options;
	options.single = true;
	options.minimizer = Options::Minimizer::Minuit2;
	
	{
		std::array<double,17> data;
		std::array<double,17> uncertainties;
		double p = rad(0);
		double a = 0;
		fill_with_sine(p,a,0,data,uncertainties);
		auto result = fit_sine(data,uncertainties,options);
		
		ASSERT_FALSE(result.valid); // Minuit2 fails in this case
		//ASSERT_APPROX(result.x[0],p,0.001);
		ASSERT_APPROX(result.x[1],a,0.001);
	}
	
	{
		std::array<double,17> data;
		std::array<double,17> uncertainties;
		double p = rad(0);
		double a = 0.02;
		fill_with_sine(p,a,0,data,uncertainties);
		auto result = fit_sine(data,uncertainties,options);
		
		ASSERT_TRUE(result.valid);
		ASSERT_TRUE(approximate_delta_periodic(result.x[0],p,0.001,AETHER_2PI));
		ASSERT_APPROX(result.x[1],a,0.001);
	}
	
	{
		std::array<double,17> data;
		std::array<double,17> uncertainties;
		double p = rad(45);
		double a = 0.02;
		fill_with_sine(p,a,0,data,uncertainties);
		auto result = fit_sine(data,uncertainties,options);
		
		ASSERT_TRUE(result.valid);
		ASSERT_APPROX(result.x[0],p,0.001);
		ASSERT_APPROX(result.x[1],a,0.001);
	}
	
	{
		std::array<double,17> data;
		std::array<double,17> uncertainties;
		double p = rad(45);
		double a = 0.04; // higher amp
		fill_with_sine(p,a,0,data,uncertainties);
		auto result = fit_sine(data,uncertainties,options);
		
		ASSERT_TRUE(result.valid);
		ASSERT_APPROX(result.x[0],p,0.001);
		ASSERT_APPROX(result.x[1],a,0.001);
	}
	
	{
		std::array<double,17> data;
		std::array<double,17> uncertainties;
		double p = rad(0);
		double a = 0.02;
		fill_with_sine(p,a,2.4,data,uncertainties);
		auto result = fit_sine(data,uncertainties,options);
		
		ASSERT_TRUE(result.valid);
		ASSERT_APPROX(result.x[0],p,0.001);
		ASSERT_APPROX(result.x[1],a,0.001);
	}
	
	{
		std::array<double,17> data;
		std::array<double,17> uncertainties;
		double p = rad(45);
		double a = 0.02;
		fill_with_sine(p,a,-2.4,data,uncertainties);
		auto result = fit_sine(data,uncertainties,options);
		
		ASSERT_TRUE(result.valid);
		ASSERT_APPROX(result.x[0],p,0.001);
		ASSERT_APPROX(result.x[1],a,0.001);
	}
	
	{
		std::array<double,17> data;
		std::array<double,17> uncertainties;
		double p = rad(-45);
		double a = 0.02;
		fill_with_sine(p,a,-2.4,data,uncertainties);
		auto result = fit_sine(data,uncertainties,options);
		
		ASSERT_TRUE(result.valid);
		ASSERT_APPROX(result.x[0],period_2pi(p),0.001);
		ASSERT_APPROX(result.x[1],a,0.001);
	}
	
	{
		std::array<double,17> data;
		std::array<double,17> uncertainties;
		double p = rad(180);
		double a = 0.02;
		fill_with_sine(p,a,-2.4,data,uncertainties);
		auto result = fit_sine(data,uncertainties,options);
		
		ASSERT_TRUE(result.valid);
		ASSERT_APPROX(result.x[0],p,0.001);
		ASSERT_APPROX(result.x[1],a,0.001);
	}
	
	{
		std::array<double,17> data;
		std::array<double,17> uncertainties;
		double p = rad(250);
		double a = 0.02;
		fill_with_sine(p,a,-2.4,data,uncertainties);
		auto result = fit_sine(data,uncertainties,options);
		
		ASSERT_TRUE(result.valid);
		ASSERT_APPROX(result.x[0],p,0.001);
		ASSERT_APPROX(result.x[1],a,0.001);
	}
}


#endif



TEST(degrees_of_freedom,test)
{

}



int main(int ,char** )
{
	RUN_ALL_TESTS();

	return EXIT_SUCCESS;
}


