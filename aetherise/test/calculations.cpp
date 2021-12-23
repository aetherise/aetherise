/**
 * \file
 *
 * \~german
 * @brief Verschiedene Berechnungen im Zusammenhang mit der Forschung.
 *
 * \~english
 * @brief Various calculations related to the research.
 *
 */


#include "test.h"

#include "../astro.h"
#include "../Theory.h"
#include "../physics.h"
#include "../mini.h"
#include "../data_reduction.h"
#include "../mathematics.h"
#include "../aetherise.h"

#include <random>
#include <sstream>
#include <functional>
#include <iomanip>

using namespace aether;





/**
 * \~german
 * Die galaktischen Koordinaten des Dipols in der kosmischen Hintergrundstrahlung
 * in das äquatoriale Koordinatensystem umwandeln, mit Fehlerfortpflanzung.
 * Dabei wird angenommen daß die Unsicherheiten von l und b nicht korreliert sind.
 *
 * \~english
 * Transform the galactic coordinates of the CMB Dipole to the equatorial coordinates
 * and propagate the uncertainties.
 * It is assumed that the uncertainties of l and b are not correlated.
 */
TEST(CMB_dipole_coordinates,galactic_to_equatorial)
{	
	Galactic gc;
	gc.l = rad(263.99);
	gc.b = rad(48.26);

	auto eq = equatorial(gc);

	Galactic gcu;
	gcu.l = rad(0.14);
	gcu.b = rad(0.03);

	real h = 0.00000001;

	{
		auto dfdl = differentiate(gc.l,[&](real x){
			return equatorial(Galactic {x,gc.b}).ra;
		},h);

		auto dfdb = differentiate(gc.b,[&](real x){
			return equatorial(Galactic {gc.l,x}).ra;
		},h);

		auto s = std::sqrt(sqr(gcu.l*dfdl) + sqr(gcu.b*dfdb));


		std::cout << "CMB dipole ra: " << rad_to_h(eq.ra) << " +/- " << rad_to_h(s) << " h\n";
	}

	{
		auto dfdl = differentiate(gc.l,[&](real x){
			return equatorial(Galactic {x,gc.b}).de;
		},h);

		auto dfdb = differentiate(gc.b,[&](real x){
			return equatorial(Galactic {gc.l,x}).de;
		},h);

		auto s = std::sqrt(sqr(gcu.l*dfdl) + sqr(gcu.b*dfdb));

		std::cout << "CMB dipole de: " << deg(eq.de) << " +/- " << deg(s) << " °\n";
	}

}



TEST(LorentzInvariance_violation,level)
{
	const double c = Lightspeed;
	const auto n = 1.00023;
	
	
	{// Theory
		auto V = CMB_dipole.v; // velocity
		auto g = sqrt_1_minus_sqr(V/c);
		auto Rd = (sqr(n)-1)/(sqr(n)+2);
		auto Rdg = Rd*g;
		auto nL = std::sqrt((1+2*Rdg)/(1-Rdg));
				
		std::cout << "(n-nL)/n = " << (n-nL)/n << "\n";
		
	}
	
	{// Result
		auto V = 326000; // velocity
		auto g = sqrt_1_minus_sqr(V/c);
		auto Rd = (sqr(n)-1)/(sqr(n)+2);
		auto Rdg = Rd*g;
		auto nL = std::sqrt((1+2*Rdg)/(1-Rdg));
				
		std::cout << "(n-nL)/n = " << (n-nL)/n << "\n";
		
	}
		
}



TEST(aether_theory,refractive_index_diff_to_amplitude_diff)
{
	const double c = Lightspeed;
	const double lat = MtWilson_Latitude;
	const auto n = 1.00023;
	
	auto V = CMB_dipole.v; // velocity
	auto g = sqrt_1_minus_sqr(V/c);
	auto Rd = (sqr(n)-1)/(sqr(n)+2);
	auto Rdg = Rd*g;
	auto nL = std::sqrt((1+2*Rdg)/(1-Rdg));
	
	std::cout << std::scientific;
	std::cout << "n-1  = " << n-1 << "\n";
	std::cout << "nL-1 = " << nL-1 << "\n";
	std::cout << "n-nL = " << n-nL << "\n";
	std::cout << std::defaultfloat;
	
	{				
		auto theory = create_theory(Options::Theory::Aether);
		auto L = Millers_Interferometer_Arm_Length;
		auto displs1 = fringe_displacements(*theory,CMB_dipole,lat,n,L,h_to_rad(17),false);
		auto displs2 = fringe_displacements(*theory,CMB_dipole,lat,nL,L,h_to_rad(17),false);
		auto amp1 = std::abs(DFT_analyze(2,displs1.begin(),displs1.end()-1));		
		auto amp2 = std::abs(DFT_analyze(2,displs2.begin(),displs2.end()-1));		
		std::cout << "amplitude difference = " << amp1-amp2 << "\n";		
	}
		
}





TEST(fringe_displacement,resting_in_ether__temperature_quadrupole)
{
	const auto c = Lightspeed;

	{
		auto n1 = 1.00023;
		auto n2 = n1 - 8e-10; // change of n, if T changes 1/1000 °C
		auto L = 32.03;
		Interferometer interferometer; // Miller's in 1925
		interferometer.light_path_1 = { {{0,L,0},n1}, {{0,-L,0},n1} }; // north
		interferometer.light_path_2 = { {{L,0,0},n2}, {{-L,0,0},n2} }; //
		interferometer.wave_length = 570e-9; // (m) acetylene lamp

		auto theory = create_theory(Options::Theory::Aether);
		Vector3 v {1,0,0}; // resting, v~0
		auto t1 = phase_delay(*theory,v,interferometer.light_path_1);
		auto t2 = phase_delay(*theory,v,interferometer.light_path_2);
		real dt = t2-t1;
		auto dl = c*dt/interferometer.wave_length;
		auto df = -2*dl;
		std::cout << "fringe displacement = " << df << "\n";
	}
}



TEST(fringe_displacement,resting_in_ether__perfect_temperature_dipole)
{
	const auto c = Lightspeed;

	{
		auto n = 1.00023;
		auto dn = 8e-8; // change of n, if T changes 1/10 °C
		auto L = 32.03;
		Interferometer interferometer; // Miller's in 1925
		interferometer.light_path_1 = { {{0,L,0},n}, {{0,-L,0},n} }; // north
		interferometer.light_path_2 = { {{L,0,0},n-dn}, {{-L,0,0},n+dn} }; //
		interferometer.wave_length = 570e-9; // (m) acetylene lamp

		auto theory = create_theory(Options::Theory::Aether);
		Vector3 v {1,0,0}; // resting, v~0
		auto t1 = phase_delay(*theory,v,interferometer.light_path_1);
		auto t2 = phase_delay(*theory,v,interferometer.light_path_2);
		auto dt = t2-t1;
		auto dl = c*dt/interferometer.wave_length;
		auto df = -2*dl;
		std::cout << "fringe displacement = " << df << "\n";
	}
}



TEST(fringe_displacement,resting_in_ether__asymmetric_temperature_dipole)
{
	const auto c = Lightspeed;

	{
		auto n = 1.00023;
		auto dn = 8e-9; // change of n, if T changes 1/100 °C
		auto L = 32.03;
		Interferometer interferometer; // Miller's in 1925
		interferometer.light_path_1 = { {{0,L,0},n}, {{0,-L,0},n} }; // north
		interferometer.light_path_2 = { {{L,0,0},n-dn}, {{-L,0,0},n+dn*1.1} }; //
		interferometer.wave_length = 570e-9; // (m) acetylene lamp

		auto theory = create_theory(Options::Theory::Aether);
		Vector3 v {1,0,0}; // resting, v~0
		auto t1 = phase_delay(*theory,v,interferometer.light_path_1);
		auto t2 = phase_delay(*theory,v,interferometer.light_path_2);
		auto dt = t2-t1;
		auto dl = c*dt/interferometer.wave_length;
		auto df = -2*dl;
		std::cout << "fringe displacement = " << df << "\n";
	}
}



TEST(fringe_displacement,unequal_arm_length)
{
	// Testet manually, even 1 m difference has no effect.
}



TEST(weather,fog)
{
	{
		auto p = barometric_formula(1700,101325,in_Kelvin(15+11.5));
		auto n1 = refractive_index_of_air(p,in_Kelvin(15),0.570,0.5,300); // clear
		auto n2 = refractive_index_of_air(p,in_Kelvin(15),0.570,1,300); // fog

		std::cout << "n1=" << n1-1 << "\n";
		std::cout << "n2=" << n2-1 << "\n";
		std::cout << "dn=" << n1-n2 << "\n";
	}
}



TEST(weather,air_pressure)
{
	{
		auto p1 = barometric_formula(1700,101000,in_Kelvin(15+11.5));
		auto p2 = barometric_formula(1700,101300,in_Kelvin(15+11.5));
		auto n1 = refractive_index_of_air(p1,in_Kelvin(15),0.570,0.5,300);
		auto n2 = refractive_index_of_air(p2,in_Kelvin(15),0.570,0.5,300);

		std::cout << "n1=" << n1-1 << "\n";
		std::cout << "n2=" << n2-1 << "\n";
		std::cout << "dn=" << n1-n2 << "\n";
	}

	{
		auto p1 = barometric_formula(1700,101300,in_Kelvin(15+11.5));
		auto p2 = barometric_formula(1700,101500,in_Kelvin(15+11.5));
		auto n1 = refractive_index_of_air(p1,in_Kelvin(15),0.570,0.5,300);
		auto n2 = refractive_index_of_air(p2,in_Kelvin(15),0.570,0.5,300);

		std::cout << "n1=" << n1-1 << "\n";
		std::cout << "n2=" << n2-1 << "\n";
		std::cout << "dn=" << n1-n2 << "\n";
	}

}



TEST(mean_index_of_refraction,Mount_Wilson)
{
	{
		auto p = barometric_formula(1700,101325,in_Kelvin(13.5+11.5));
		auto n = refractive_index_of_air(p,in_Kelvin(13.5),0.570,0.5,305);
		std::cout << "mean p = " << p << "\n";
		std::cout << "mean n-1 = " << n-1 << "\n";
	}
}



TEST(change_of_signal_while_measuring,error)
{	
	const double lat = MtWilson_Latitude;
	const auto n = 1.00023;
	const double L = Millers_Interferometer_Arm_Length;
	const int N = 20; // number of turns
	//const double sidereal_time = 0; // h
	const double dt = (15./60.)/N; // 15 mins, N turns, time of 1 turn in h
	
	double max_dA = 0;
	double max_dphi = 0;
	double max_ui=0;
	
	for(int k=0;k<24;k++) 
	{			
		double sidereal_time = k; // h
		std::vector<std::array<double,17>> turns;
		
		auto theory = create_theory(Options::Theory::Aether);
		for (int i=0;i<N;i++) {
			auto theta = h_to_rad(sidereal_time+i*dt);
			turns.push_back(fringe_displacements(*theory,CMB_dipole,lat,n,L,theta,false));			
		}
		
		// mean displacements
		std::array<double,17> means {};
		for (const auto& turn : turns) {
			add_array(means,turn);
		}		
		for (auto& mean : means) {
			mean /= turns.size();
		}
		
		// mean signal
		auto z = DFT_analyze(2,means.begin(),means.end()-1);		
		auto A = std::abs(z);
		auto phi = std::arg(z);
		
		// signal at mean observation time
		auto theta = h_to_rad(sidereal_time+(N-1)*dt/2.);
		auto mot_displs = fringe_displacements(*theory,CMB_dipole,lat,n,L,theta,false);		
		auto mot_z = DFT_analyze(2,mot_displs.begin(),mot_displs.end()-1);		
		auto mot_A = std::abs(mot_z);
		auto mot_phi = std::arg(mot_z);
		
		//std::cout << "(A, phi)_mean = (" << A << ", " << deg(phi) << "°)\n";		
		//std::cout << "(A, phi)_mot  = (" << mot_A << ", " << deg(mot_phi) << "°)\n";		
		//std::cout << "(dA, dphi)    = (" << A-mot_A << ", " << deg(phi-mot_phi) << "°)\n";		
		//std::cout << "k=" << k << "\n";
		if (std::abs(A-mot_A) > std::abs(max_dA))
			max_dA = A-mot_A;
		if (std::abs(phi-mot_phi) > std::abs(max_dphi))
			max_dphi = phi-mot_phi;
				
		
		// Uncertainty
				
		std::array<double,17> uis;
		for (size_t i=0;i<17;i++) {
			std::vector<double> qis;
			for (const auto& turn : turns) {
				qis.push_back(turn.at(i));
			}
			auto s = sample_standard_deviation(qis.begin(),qis.end(),means.at(i));
			uis.at(i) = s/std::sqrt(turns.size());			 
		}
		
		for (const auto& ui : uis)
			max_ui = std::max(max_ui,ui);
		
	}
	
	std::cout << "(dA, dphi)_max = (" << max_dA << ", " << rad_to_h(max_dphi) << " h)\n";		
	std::cout << "max u = " << max_ui << "\n";		
}




TEST(shankland1955,sidereal_times)
{
	{
		// To test Fig. 4(A)
		
		// Cleveland, Ohio
		//auto lat = rad(41.504);
		auto lon = rad(-81.608);
		auto tz = -5.; // timezone 
		
		Calendar cal = Calendar {1927,8,30 + (0-tz)/24.}; // be carefull with tz: overflow or neg. values not allowed
		auto theta = sidereal_time(cal,lon);
		std::cout << "sidereal time on 1927-8-30 at 0:00 in Cleveland: " << h_to_time(rad_to_h(theta)) << "\n";		
	}
	
	{
		// Cleveland, Ohio
		//auto lat = rad(41.504);
		auto lon = rad(-81.608);
		auto tz = -5.; // timezone 
		
		Calendar cal = Calendar {1924,7,8 + (12-tz)/24.}; 
		auto theta = sidereal_time(cal,lon);
		std::cout << "sidereal time on 1924-7-8 at 12:00 in Cleveland: " << h_to_time(rad_to_h(theta)) << "\n";		
	}
}



TEST(shankland1955,index_of_refraction)
{
	{
		double p = 101325; // Pa
		double T = in_Kelvin(18); // °C -> K
		double lambda = Millers_Interferometer_Wave_Length*1e+6; // µm 
		double h = 0.5; // 50%
		double xc = 305.8; // ppm 
		double n = refractive_index_of_air(p,T,lambda,h,xc);
		//std::cout << std::setprecision(16);
		std::cout << "Index of Refraction in Cleveland on 1927-08-30 at 18°C: " << n << "\n";
	}
}



TEST(Pease1930,sidereal_times)
{
	{
		// 
		
		// Pasadena
		//auto lat = 
		auto lon = MtWilson_Longitude; // near enough
		auto tz = -8.; // timezone 
		
		Calendar cal = Calendar {1928,1,6 + (12-tz)/24.}; // be carefull with tz: overflow or neg. values not allowed
		auto theta = sidereal_time(cal,lon);
		std::cout << "sidereal time on 1928-1-6 at 12:00 in Pasadena: " << h_to_time(rad_to_h(theta)) << "\n";		
	}
	
	{
		// 
		
		// Pasadena
		//auto lat = 
		auto lon = MtWilson_Longitude; // near enough
		auto tz = -8.; // timezone 
		
		Calendar cal = Calendar {1928,1,13 + (22-tz)/24.}; // be carefull with tz: overflow or neg. values not allowed
		auto theta = sidereal_time(cal,lon);
		std::cout << "sidereal time on 1928-1-13 at 22:00 in Pasadena: " << h_to_time(rad_to_h(theta)) << "\n";		
	}
	
}



TEST(propagation_of_uncertainty,pre_processed_samples)
{
	const auto seed = create_random_seed();
	const int N = 100;
	
	std::mt19937 engine(seed);
	std::normal_distribution<double> dist(0,0.4);

	std::vector<double> samples1;		
	for (int i=0;i<N;i++) {
		samples1.push_back(dist(engine));			
	}
	std::vector<double> samples2;
	for (int i=0;i<N;i++) {
		samples2.push_back(dist(engine));			
	}
	
	{		
		std::cout << "mean of estimates: error propagation\n";
		auto mean1 = mean_value(samples1.begin(),samples1.end());
		auto sigma1 = sample_standard_deviation(samples1.begin(),samples1.end(),mean1);
		auto u1 = sigma1/std::sqrt(samples1.size());
		
		auto mean2 = mean_value(samples2.begin(),samples2.end());
		auto sigma2 = sample_standard_deviation(samples2.begin(),samples2.end(),mean2);
		auto u2 = sigma2/std::sqrt(samples2.size());

		auto x = (mean1+mean2)*0.5;
		auto u = 0.5*std::sqrt(sqr(u1)+sqr(u2));
		std::cout << "x=" << x << ", u=" << u << "\n";
	}
	
	{				
		std::cout << "sample mean: uncertainty\n";
		std::vector<double> samples;
		for (size_t i=0;i<samples1.size();i++) {
			samples.push_back(0.5*(samples1.at(i)+samples2.at(i)));
		}
		auto mean = mean_value(samples.begin(),samples.end());
		auto sigma = sample_standard_deviation(samples.begin(),samples.end(),mean);
		auto u = sigma/std::sqrt(samples.size());
				
		std::cout << "x=" << mean << ", u=" << u << "\n";
	}
}



TEST(phase_with_variance,aggregate_comparison)
{
	const auto seed = create_random_seed();
	
	{		
		std::cout << "double period\n";
		std::mt19937 engine(seed);
		std::normal_distribution<double> dist(0,0.4);
		
		std::array<double,17> samples;
		double sumA = 0;
		double sumphi = 0;
		double sumZA = 0;
		double sumZphi = 0;
		int N = 100;
		std::complex<double> z {};
		DFTGoertzel dft({2},16);
		for (int i=0;i<N;i++) {
			double A = 1;	
			double phi = 0+dist(engine);
			set_sine(A,phi,0,8,samples);
			add_sine(1,2,0,16,samples);
			dft.analyze(samples.begin(),samples.end()-1);
			auto Z = DFT_analyze(2,samples.begin(),samples.end()-1);
			sumZA += std::abs(Z);
			sumZphi += std::arg(Z)+AETHER_PI_2;
			
			z += std::polar(A,phi);
			sumA += A;
			sumphi += phi;
		}
		auto Z = dft.result().at(2);
		
		// superimposed				
		z /= double(N);				
		std::cout << "A=" << std::abs(z) << ", phi=" << std::arg(z) << "\n";
		std::cout << "A=" << std::abs(Z) << ", phi=" << std::arg(Z)+AETHER_PI_2 << " (DFT)\n";
		// not superimposed
		std::cout << "A=" << sumA/N      << ", phi=" << sumphi/N << " (mean)\n";
		std::cout << "A=" << sumZA/N     << ", phi=" << sumZphi/N << " (DFT mean)\n";
		
	}
	
	{		
		std::cout << "2x half double period\n";
		std::mt19937 engine(seed);
		std::normal_distribution<double> dist(0,0.4);
		
		std::array<double,17> samples;
		double sumA = 0;
		double sumphi = 0;
		double sumZA = 0;
		double sumZphi = 0;
		int N = 100;
		std::complex<double> z {};
		DFTGoertzel dft({1},8);
		for (int i=0;i<N;i++) {
			double A = 1;	
			double phi = 0+dist(engine);
			set_sine(A,phi,0,8,samples);
			add_sine(1,2,0,16,samples); // without this, sum of (A,phi) works
			dft.analyze(samples.begin(),samples.end()-1-8);
			dft.analyze(samples.begin()+8,samples.end()-1);
			auto Z = DFT_analyze(1,samples.begin(),samples.end()-1-8);
			sumZA += std::abs(Z);
			sumZphi += std::arg(Z)+AETHER_PI_2;
			Z = DFT_analyze(1,samples.begin()+8,samples.end()-1);
			sumZA += std::abs(Z);
			sumZphi += std::arg(Z)+AETHER_PI_2;
			
			z += std::polar(A,phi);
			sumA += A;
			sumphi += phi;
		}
		auto Z = dft.result().at(1);
		
		
		// superimposed				
		z /= double(N);				
		std::cout << "A=" << std::abs(z) << ", phi=" << std::arg(z) << "\n";
		std::cout << "A=" << std::abs(Z) << ", phi=" << std::arg(Z)+AETHER_PI_2 << " (DFT)\n";
		// not superimposed
		std::cout << "A=" << sumA/N      << ", phi=" << sumphi/N << " (mean)\n";
		std::cout << "A=" << sumZA/(N*2) << ", phi=" << sumZphi/(N*2) << " (DFT mean)\n";
		
	}
}



TEST(sidereal_times,validation)
{
	{
		// Apr-110
		// TODO why not +8 min? Millers sidereal time incorrect?
		Calendar cal {1925,3,31 + time_to_h({00,38})/24. /*+ (8/60.)/24.*/}; 
		auto theta = sidereal_time(cal,0);
		std::cout << "Apr-110: noted 13:11, calculated " << h_to_time(rad_to_h(theta)) << "\n";
		
	}
	
	{
		// Aug-60
		Calendar cal {1925,8,5 + time_to_h({13,10})/24. + (8/60.)/24.};
		auto theta = sidereal_time(cal,0);
		std::cout << "Aug-60: noted 10:10, calculated " << h_to_time(rad_to_h(theta)) << "\n";
		
	}
	
	{
		// Sep-49
		Calendar cal {1925,9,17 + time_to_h({19,30})/24. + (8/60.)/24.};
		auto theta = sidereal_time(cal,0);
		std::cout << "Sep-49: noted 19:23, calculated " << h_to_time(rad_to_h(theta)) << "\n";
		
	}
	
	{
		// Feb-43
		Calendar cal {1926,2,7 + time_to_h({17,27})/24. + (8/60.)/24.};
		auto theta = sidereal_time(cal,0);
		std::cout << "Feb-43: noted 2:43, calculated " << h_to_time(rad_to_h(theta)) << "\n";
		
	}
}





int main(int,char**)
{
	RUN_ALL_TESTS();

	return EXIT_SUCCESS;
}
