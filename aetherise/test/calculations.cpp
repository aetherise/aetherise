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
		auto dfdl = derive(gc.l,[&](real x){
			return equatorial(Galactic {x,gc.b}).ra;
		},h);

		auto dfdb = derive(gc.b,[&](real x){
			return equatorial(Galactic {gc.l,x}).ra;
		},h);

		auto s = std::sqrt(sqr(gcu.l*dfdl) + sqr(gcu.b*dfdb));


		std::cout << "CMB dipole ra: " << rad_to_h(eq.ra) << " +/- " << rad_to_h(s) << " h\n";
	}

	{
		auto dfdl = derive(gc.l,[&](real x){
			return equatorial(Galactic {x,gc.b}).de;
		},h);

		auto dfdb = derive(gc.b,[&](real x){
			return equatorial(Galactic {gc.l,x}).de;
		},h);

		auto s = std::sqrt(sqr(gcu.l*dfdl) + sqr(gcu.b*dfdb));

		std::cout << "CMB dipole de: " << deg(eq.de) << " +/- " << deg(s) << " °\n";
	}

}



TEST(aether_theory,refractive_index_diff_to_amplitude_diff)
{
	{
		auto theory = create_theory(Options::Theory::Aether);
		auto displs1 = fringe_displacements(*theory,CMB_dipole,1.00023,14,false);
		auto displs2 = fringe_displacements(*theory,CMB_dipole,1.00022,14,false);
		auto amp1 = max_abs_value(displs1);
		auto amp2 = max_abs_value(displs2);
		std::cout << "amplitude difference = " << amp1-amp2 << "\n";
		std::cout << "amplitude ratio = " << amp1/amp2 << "\n";
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




#ifdef AETHER_MINUIT


template <size_t N>
class FunctionContext : public MinimizeContext
{
	std::function<double(const std::vector<double>&)> f;
public:
	FunctionContext(const std::function<double(const std::vector<double>&)>& f)
		:f{f} {
	}

	double operator() (const std::vector<double>& params) override
	{
		return f(params);
	}


	double delta_chi_squared() const {
		return 1.0;
	}

	bool fixed_ad() const {
		return false;
	}
};




TEST(chi_squared_test,pre_processed_data)
{
	{
		constexpr int N = 17;
		const int n = 100;

		auto fx = [](double x){
			return 0.5*x*x + 2*x - 1;
		};
		auto fx2 = [](double x){
			return 0.5*std::pow(x,2) + 2*x - 1.5;
		};


		std::mt19937 engine(create_random_seed());
		std::normal_distribution<double> dist(0,1);

		double jX = 0;
		int jn = 0;
		for (int j=0;j<100;j++) {
			std::vector<std::array<double,N>> as;
			std::vector<std::array<double,N>> us;
			for (int ia=0;ia<17;ia++) {
				std::array<double,N> a;
				std::array<double,N> u;

				for (int i=0;i<N;i++) {
					auto y = ia==0 && i<N/2 ? fx2(i) : fx(i);

					std::vector<double> m;
					for (int k=0;k<n;k++) {
						m.push_back(y + dist(engine));
					}

					// estimate
					double s = 0;
					for (auto&& mi : m) {
						s += mi;
					}
					s/=n;
					a.at(i) = s;

					// uncertainty
					s = 0;
					for (auto&& mi : m) {
						s += sqr(mi-a.at(i));
					}
					s/=n*(n-1);
					u.at(i) = std::sqrt(s);
				}
				as.push_back(std::move(a));
				us.push_back(std::move(u));
			}


			FunctionContext<N> context([&](const std::vector<double>& params){
				std::array<double,N> t;
				for (int i=0;i<N;i++) {
					t.at(i) = params.at(0)*i*i + params.at(1)*i + params.at(2);
				}



				std::array<double,N> a {};
				std::array<double,N> u {};
/*

				const auto ni = 3;
				for (int i=0;i<ni;i++) {
					add_array(a,as.at(0));
					add_sqr_array(u,us.at(i));
				}
				for (auto&& ai : a) {
					ai/=ni;
				}
				for (auto&& ui : u) {
					ui = std::sqrt(ui)/ni;
				}*/

				double chi=0;
/*
				add_array(a,as.at(0));
				add_array(a,as.at(1));
				add_sqr_array(u,us.at(0));
				add_sqr_array(u,us.at(1));
				for (auto&& ai : a)
					ai/=2;
				for (auto&& ui : u)
					ui = std::sqrt(ui)/2;
				chi += chi_squared_test(a,u,t);

				a = {};
				u = {};
				add_array(a,as.at(1));
				add_array(a,as.at(2));
				add_sqr_array(u,us.at(1));
				add_sqr_array(u,us.at(2));
				for (auto&& ai : a)
					ai/=2;
				for (auto&& ui : u)
					ui = std::sqrt(ui)/2;
				chi += chi_squared_test(a,u,t);

				a = {};
				u = {};
				add_array(a,as.at(2));
				add_array(a,as.at(0));
				add_sqr_array(u,us.at(2));
				add_sqr_array(u,us.at(0));
				for (auto&& ai : a)
					ai/=2;
				for (auto&& ui : u)
					ui = std::sqrt(ui)/2;
				chi += chi_squared_test(a,u,t);
	*/


	/*			t={};

				add_array(a,as.at(0));
				sub_array(a,as.at(1));
				add_sqr_array(u,us.at(0));
				add_sqr_array(u,us.at(1));
				for (auto&& ui : u)
					ui = std::sqrt(ui);
				chi += chi_squared_test(a,u,t);

				a = {};
				u = {};
				add_array(a,as.at(2));
				sub_array(a,as.at(3));
				add_sqr_array(u,us.at(2));
				add_sqr_array(u,us.at(3));
				for (auto&& ui : u)
					ui = std::sqrt(ui);
				chi += chi_squared_test(a,u,t);

				a = {};
				u = {};
				add_array(a,as.at(4));
				sub_array(a,as.at(5));
				add_sqr_array(u,us.at(4));
				add_sqr_array(u,us.at(5));
				for (auto&& ui : u)
					ui = std::sqrt(ui);
				chi += chi_squared_test(a,u,t);

*/


			/*	chi=0;
				for (size_t i=0;i<9;i++)
					chi += chi_squared_test(as.at(i),us.at(i),t,9);*/


				reduce_to_single_period(t);
				chi=0;
				for (size_t i=0;i<9;i++) {
					ReducedData rd;
					rd.displacements = as.at(i);
					rd.uncertainties = us.at(i);
					reduce_to_single_period(rd);
					chi += chi_squared_test(rd.displacements,rd.uncertainties,t,9);
				}

				chi=0;
				auto t2 = t;
				add_array(t2,t);
				for (size_t i=0;i<9;i++) {
					a = as.at(i);
					add_array(a,as.at((i+1)%9));
					u={};
					add_sqr_array(u,us.at(i));
					add_sqr_array(u,us.at((i+1)%9));
					for (auto&& ui : u)
						ui = std::sqrt(ui);
					chi += chi_squared_test(a,u,t2,9);
				}

				return chi;
			});

			std::ostringstream ss;
			auto result = minimize_locally_Minuit2({0.5,2,-1},context,ss);
			if (!result.valid) {
				std::cerr << ss.str();
				continue;
			}

			jX += result.y;
			jn++;
		}

		std::cout << "expected dof = " << 9*(9)-3 << "\n";
		std::cout << "jn=" << jn << "\n";
		std::cout << "X²=" << jX/jn << "\n";
	}
}




#endif // Minuit dependancy



TEST(azimuth_measurement,test_for_normality)
{
	// https://de.wikipedia.org/wiki/Konfidenzintervall_f%C3%BCr_die_Erfolgswahrscheinlichkeit_der_Binomialverteilung
	{
		int n = 1224;
		int k = 1125;
		//double p = 0.93;

		double c = 1.96; // 95% confidence

		auto k_= k+c*c/2;
		auto n_ = n + c*c;
		auto p_ = k_/n_;

		auto u_ = c*std::sqrt(p_*(1-p_)/n_);

		std::cout << "p_ = " << p_ << " +/- " << u_ << "\n";
		std::cout << "real nnd: " << 100*(0.95-p_) << "% +/- " << u_*100 << "%\n";
	}
}



TEST(azimuth_measurement,interval_test)
{
	{
		std::mt19937 engine(create_random_seed());
		std::normal_distribution<double> dist(0,0.02);

		int m = 0;
		for (int j=0;j<100;j++) {
			std::vector<double> x;

			const int n = 1000;
			int nnd = 0;
			for (int i=0;i<n;i++) {
				x.clear();
				if (i<n-n*0.02) { // 2% not normal
					// normal
					for (int k=0;k<20;k++) {
						x.push_back(dist(engine));
					}
				}
				else {
					// double normal: systematic error change
					for (int k=0;k<10;k++) {
						x.push_back(dist(engine));
					}
					for (int k=0;k<10;k++) {
						x.push_back(dist(engine)+0.2);
					}
				}
				std::sort(x.begin(),x.end());
				auto A = test_for_normality(x);
				if (A > ADTestDAgostinoQuantiles.q_5)
					nnd++;
			}

			ASSERT_APPROX(nnd,70,0.5);
			auto k = n-nnd;

			double c = 1.96; // 95% confidence
			auto k_= k+c*c/2;
			auto n_ = n + c*c;
			auto p_ = k_/n_;
			auto u_ = c*std::sqrt(p_*(1-p_)/n_);

			const double p = 0.93;
			if (p>=p_-u_ && p<=p_+u_)
				m++;
		}

		std::cout << "inside interval: " << m << "\n";
	}
}



TEST(change_of_signal_while_measuring,uncertainty)
{
	// Variance: s² = sum((xi-x)²)
	// so, with d=0.002/20=0.0001 and x=0 and xi=i
	// s² = sum(d²*(i-n/2)²)
	double d = 0.0001;
	double n = 20; // turns
	auto s = std::sqrt( d*d * n*n * (1./3 - 1./2 + 1./4));
	auto u = s/std::sqrt(20);
	
	std::cout << "u=" << u << "\n";
}



template<size_t N>
void fill_with_sine(double p, double a,double c,std::array<double,N>& data)
{
	for (size_t i=0;i<N;i++) {
		data.at(i) = a*std::sin(i*AETHER_PI/4-p) + c;		
	}
}





TEST(degrees_of_freedom,phase_amplitude)
{	
	const double sigma_p = 0.3*4;
	const double sigma_a = 0.005*4;
	std::mt19937 engine(create_random_seed());
	std::normal_distribution<double> distp(0,sigma_p);
	std::normal_distribution<double> dista(0,sigma_a);
	std::cout << "sigma p = " << sigma_p << "\n";
	std::cout << "sigma a = " << sigma_a << "\n";
	Options options;
	options.single = true;
	
	const double a = 0.02;		
	const int M = 32;
	const int N = 40;
	const int TURNS = 20;
	std::cout << "a = " << a << "\n";	
	
	double chi = 0;
	double min_chi = 0;
	for (int x=0;x<M;x++) {
		const double p = AETHER_2PI/M*x;
		
		std::array<double,8> theory;
		fill_with_sine(p,a,0,theory);
		//std::cout << "p = " << p << "\n";
			
		for (int j=0;j<N;j++) {
			std::vector<std::array<double,8>> turns;	
			for (int k=0;k<TURNS;k++) {
				std::array<double,8> data;		
				fill_with_sine(p+distp(engine),a+dista(engine),0,data);	
				turns.push_back(std::move(data));
			}
			
			std::array<double,8> data;
			std::array<double,8> uncertainties;	
			for (int i=0;i<8;i++) {
				std::vector<double> s;
				for (auto& turn : turns)
					s.push_back(turn.at(i));
				auto mean = mean_value(s.begin(),s.end());
				data.at(i) = mean;
				uncertainties.at(i) = sample_standard_deviation(s.begin(),s.end(),mean);
			}
			for (int i=0;i<8;i++)
				uncertainties.at(i) /= std::sqrt(turns.size());
			
			/*std::cout << "data: ";
			output_separated(std::cout,data.begin(),data.end(),", ");
			std::cout << "\n";*/
			/*std::cout << "uncertainties: ";
			output_separated(std::cout,uncertainties.begin(),uncertainties.end(),", ");
			std::cout << "\n";*/
				
			chi += chi_squared_test(data,uncertainties,theory);
			
			std::array<double,17> data17;
			std::array<double,17> uncertainties17;
			for (size_t i=0;i<8;i++) {
				data17.at(i) = data.at(i);
				data17.at(i+8) = data.at(i);						
				uncertainties17.at(i) = uncertainties.at(i);
				uncertainties17.at(i+8) = uncertainties.at(i);
			}
			data17.at(16) = data17.at(0);
			uncertainties17.at(16) = uncertainties17.at(0);
			
			auto local_min = fit_sine(data17,uncertainties17,options);
			
			min_chi += local_min.y;
		}	
	}
	std::cout << "X² = " << chi/(N*M) << "\n";	
	std::cout << "min X² = " << min_chi/(N*M) << "\n";	
}



TEST(normal_distribution,test)
{	
	std::mt19937 engine(create_random_seed());	
	std::normal_distribution<double> dist(3.1,0.01);
	
	const int n = 100;
	std::vector<double> values;
	for (int i=0;i<n;i++) {
		values.push_back(dist(engine));
	}
	
	auto mean = mean_value(values.begin(),values.end());
	auto ssd = sample_standard_deviation(values.begin(),values.end(),mean);
	
	ASSERT_APPROX(mean,3.1,0.2);
	ASSERT_APPROX(ssd,0.01,0.2);
}





TEST(shankland1955,sidereal_times)
{
	{
		// To test Fig. 4(A)
		
		// Cleveland, Ohio
		//auto lat = 41.497118;
		auto lon = -81.679100;
		auto tz = -5.; // timezone -5
		
		Calendar cal = Calendar {1927,8,30 + (0-tz)/24.}; 
		auto theta = sidereal_time(cal,lon);
		std::cout << "sidereal time in 1927-8-30 at 0:00 in Cleveland: " << h_to_time(rad_to_h(theta)) << "\n";		
	}
	
	{
		// Cleveland, Ohio
		//auto lat = 41.497118;
		auto lon = -81.679100;		
		auto tz = -5.; // timezone -5
		
		Calendar cal = Calendar {1924,7,8 + (12-tz)/24.}; 
		auto theta = sidereal_time(cal,lon);
		std::cout << "sidereal time in 1924-7-8 at 12:00 in Cleveland: " << h_to_time(rad_to_h(theta)) << "\n";		
	}
}



int main(int,char**)
{
	RUN_ALL_TESTS();

	return EXIT_SUCCESS;
}
