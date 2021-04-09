#include "test.h"

#include "../mini.h"
#include "../Theory.h"
#include "../utils.h"

#include <random>

using namespace aether;


const double Ri = 1.00023; // refractive index

TEST(minimize_theory,without_error)
{
	auto theory = create_theory(Options::Theory::Aether);
	const double f = 100000;
	const int max_iter = 20;

	std::array<double,3> x0 {250000.0/f,h_to_rad(16.0),rad(7.0)};
	std::array<double,17> uncertainties;
	std::fill(uncertainties.begin(),uncertainties.end(),0.01);

	auto local_min = minimize_locally_a(x0,[&](const std::array<double,3>& x){
		double chi = 0;
		for (int i=0;i<24;i++) {
			auto sidereal_time = h_to_rad(i);
			auto displs = fringe_displacements(*theory,CMB_dipole,Ri,sidereal_time,false);
			TheoryParameters params {x[0]*f,x[1],x[2]};
			auto xdispls = fringe_displacements(*theory,params,Ri,sidereal_time,false);
			chi += chi_squared_test(displs,uncertainties,xdispls);
		}
		return chi;
	},0.000001,max_iter,0.000001);

	ASSERT_TRUE(local_min.i<max_iter);

	TheoryParameters min_params {local_min.x[0]*f,local_min.x[1],local_min.x[2]};
	std::cout << "v=" << min_params.v << "\n";
	std::cout << "ra=" << rad_to_h(min_params.a) << "\n";
	std::cout << "de=" << deg(min_params.d) << "\n";

	ASSERT_APPROX(min_params.v,CMB_dipole.v,0.01);
	ASSERT_APPROX(min_params.a,CMB_dipole.a,0.01);
	ASSERT_APPROX(min_params.d,CMB_dipole.d,0.01);
}



TEST(minimize_theory,random_error)
{
	auto theory = create_theory(Options::Theory::Aether);
	const double f = 100000;
	const int max_iter = 20;


	auto seed = create_random_seed();
	std::mt19937 engine(seed);
	std::uniform_real_distribution<double> dist(-1,1);


	std::array<double,3> x0 {250000.0/f,h_to_rad(16.0),rad(7.0)};
	std::array<double,17> uncertainties;
	std::fill(uncertainties.begin(),uncertainties.end(),0.01);

	auto local_min = minimize_locally_a(x0,[&](const std::array<double,3>& x){
		engine.seed(seed);
		double chi = 0;
		for (int i=0;i<24;i++) {
			auto sidereal_time = h_to_rad(i);
			auto displs = fringe_displacements(*theory,CMB_dipole,Ri,sidereal_time,false);
			for (auto& d: displs) {
				d += dist(engine)/200;
			}
			TheoryParameters params {x[0]*f,x[1],x[2]};
			auto xdispls = fringe_displacements(*theory,params,Ri,sidereal_time,false);

			chi += chi_squared_test(displs,uncertainties,xdispls);
		}
		return chi;
	},0.000001,max_iter,0.000001);

	std::cout << "i=" << local_min.i << "\n";
	//ASSERT_TRUE(local_min.i<max_iter);

	TheoryParameters min_params {local_min.x[0]*f,local_min.x[1],local_min.x[2]};
	std::cout << "v=" << min_params.v << "\n";
	std::cout << "ra=" << rad_to_h(min_params.a) << "\n";
	std::cout << "de=" << deg(min_params.d) << "\n";


	//ASSERT_APPROX(min_params.v,CMB_dipole.v,0.1);
	//ASSERT_APPROX(min_params.a,CMB_dipole.a,0.1);
	//ASSERT_APPROX(min_params.d,CMB_dipole.d,0.1);
}



TEST(minimize_theory,systematic_error)
{
	auto theory = create_theory(Options::Theory::Aether);
	const double f = 100000;
	const int max_iter = 20;


	auto seed = create_random_seed();
	std::mt19937 engine(seed);
	std::uniform_real_distribution<double> dist(-1.0,1.0);


	std::array<double,3> x0 {250000.0/f,h_to_rad(16.0),rad(7.0)}; // should be searched for
	std::array<double,17> uncertainties;
	std::fill(uncertainties.begin(),uncertainties.end(),0.01);

	auto local_min = minimize_locally_a(x0,[&](const std::array<double,3>& x){
		engine.seed(seed);
		double chi = 0;
		for (int i=0;i<24;i++) {
			auto sidereal_time = h_to_rad(i);
			auto displs = fringe_displacements(*theory,CMB_dipole,Ri,sidereal_time,false);

			auto phase = dist(engine)/2;
			auto amplitude = 0.01*(2+dist(engine));
			int t=0;
			for (auto& d: displs) {
				d += amplitude*sin(AETHER_PI/8.0*t-phase);
				t++;
			}
			TheoryParameters params {x[0]*f,x[1],x[2]};
			auto xdispls = fringe_displacements(*theory,params,Ri,sidereal_time,false);

			chi += chi_squared_test(displs,uncertainties,xdispls);
		}
		return chi;
	},0.000001,max_iter,0.000001);

	std::cout << "i=" << local_min.i << "\n";
	//ASSERT_TRUE(local_min.i<max_iter);

	TheoryParameters min_params {local_min.x[0]*f,local_min.x[1],local_min.x[2]};
	std::cout << "v=" << min_params.v << "\n";
	std::cout << "ra=" << rad_to_h(min_params.a) << "\n";
	std::cout << "de=" << deg(min_params.d) << "\n";


	//ASSERT_APPROX(min_params.v,CMB_dipole.v,0.1);
	//ASSERT_APPROX(min_params.a,CMB_dipole.a,0.1);
	//ASSERT_APPROX(min_params.de,CMB_dipole.de,0.1);
}




TEST(minimize_theory,random_and_systematic_error)
{
	auto theory = create_theory(Options::Theory::Aether);
	const double f = 100000;
	const int max_iter = 20;


	auto seed = create_random_seed();
	std::mt19937 engine(seed);
	std::uniform_real_distribution<double> dist(-1.0,1.0);


	std::array<double,3> x0 {250000.0/f,h_to_rad(16.0),rad(7.0)}; // should be searched for
	std::array<double,17> uncertainties;
	std::fill(uncertainties.begin(),uncertainties.end(),0.01);

	auto local_min = minimize_locally_a(x0,[&](const std::array<double,3>& x){
		engine.seed(seed);
		double chi = 0;
		for (int i=0;i<24;i++) {
			auto sidereal_time = h_to_rad(i);
			auto displs = fringe_displacements(*theory,CMB_dipole,Ri,sidereal_time,false);

			auto phase = dist(engine)/2;
			auto amplitude = 0.01*(2+dist(engine));
			int t=0;
			for (auto& d: displs) {
				d += amplitude*sin(AETHER_PI/8.0*t-phase); // systematic
				d += dist(engine)/200; // random
				t++;
			}
			TheoryParameters params {x[0]*f,x[1],x[2]};
			auto xdispls = fringe_displacements(*theory,params,Ri,sidereal_time,false);


			chi += chi_squared_test(displs,uncertainties,xdispls);
		}
		return chi;
	},0.000001,max_iter,0.000001);

	std::cout << "i=" << local_min.i << "\n";
	//ASSERT_TRUE(local_min.i<max_iter);

	TheoryParameters min_params {local_min.x[0]*f,local_min.x[1],local_min.x[2]};
	std::cout << "v=" << min_params.v << "\n";
	std::cout << "ra=" << rad_to_h(min_params.a) << "\n";
	std::cout << "de=" << deg(min_params.d) << "\n";


	//ASSERT_APPROX(min_params.v,CMB_dipole.v,0.1);
	//ASSERT_APPROX(min_params.a,CMB_dipole.a,0.1);
	//ASSERT_APPROX(min_params.de,CMB_dipole.de,0.1);
}




int main(int ,char** )
{
	RUN_ALL_TESTS();

	return EXIT_SUCCESS;
}


