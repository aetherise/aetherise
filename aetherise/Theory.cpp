#include "Theory.h"
#include "mathematics.h"
#include "astro.h"
#include "utils.h"



namespace aether {





//--------------
// Theories
//--------------

real Classic::phase_delay(const Vector3& v,const Lightpath& lp) const
{
	const real c = Lightspeed;

	auto V = length(v); // velocity

	auto x = dot_product(v,lp.p)/V;
	auto px = v/V * x;
	auto py = lp.p-px;
	auto y = length(py); // losing sign

	auto a = x*V/sub_sqr(c/lp.n,V);
	auto b = (x*x+y*y)/sub_sqr(c/lp.n,V);
	return a + std::sqrt( sqr(a)+b );
}




real Ether::phase_delay(const Vector3& v,const Lightpath& lp) const
{
	const real c = Lightspeed;

	auto V = length(v); // velocity
	auto g = sqrt_1_minus_sqr(V/c);

	auto x = dot_product(v,lp.p)/V;
	auto px = v/V * x;
	auto py = lp.p-px;
	auto y = length(py); // losing sign


	// anisotropic index of refraction
	auto Rd = (sqr(lp.n)-1)/(sqr(lp.n)+2);
	auto Rdg = Rd*g;
	auto ng = std::sqrt((1+2*Rdg)/(1-Rdg));

	auto L = length(lp.p);
	auto ux_ = x/L*(c/ng);
	auto uy_ = y/L*(c/lp.n);

	// velocity-addition formula		
	auto ux = (ux_+V) / (1 + (ux_/c)*(V/c));
	auto uy = (uy_*g) / (1 + (ux_/c)*(V/c));
	auto u = std::sqrt(sqr(ux-V) + sqr(uy));

	auto Lg = std::sqrt(sqr(x*g) + sqr(y));

	return Lg/u;
}



real Relativity::phase_delay(const Vector3& ,const Lightpath& lp) const
{
	const real c = Lightspeed;

	auto L = length(lp.p);
	return L/(c/lp.n);
}




//-------------------
// Functions
//-------------------



real phase_time(const Theory& theory,const Vector3& v,const std::vector<Lightpath>& light_path)
{
	real t = 0;
	for (const Lightpath& lp : light_path) {
		t += theory.phase_delay(v,lp);
	}

	return t;
}



std::unique_ptr<Theory> create_theory(Options::Theory options_theory)
{
	std::unique_ptr<Theory> theory;

	switch(options_theory) {
	case Options::Theory::Classic:
		theory = make_unique<Classic>();
		break;	
	case Options::Theory::Ether:
		theory = make_unique<Ether>();
		break;
	case Options::Theory::Relativity:
		theory = make_unique<Relativity>();
		break;	
	default:
		throw std::runtime_error("unknown theory");
	}

	return theory;
}




std::array<double,17>
fringe_displacements(const Theory& theory,const TheoryParameters& params,double n,
					 double sidereal_time,bool invert)
{
	const real c = Lightspeed;

	Equatorial apex {params.a, params.d};
	Horizontal hz = horizontal(resting(apex,sidereal_time),rad(34.225)); // latitude Mt. Wilson

	Vector3 v {0,params.v,0}; // (m/s) velocity vector, pointing north
	v = rotate_x(v,hz.h); // rotate to be parallel to the given apex
	v = rotate_z(v,-north_azimuth(hz.a));

	real L = 32.03; // (m) virtual arm length of the Michelson interferometer	
	Interferometer interferometer; // Miller's in 1925
	interferometer.light_path_1 = { {{0,L,0},n}, {{0,-L,0},n} }; // north
	interferometer.light_path_2 = { {{L,0,0},n}, {{-L,0,0},n} }; //
	interferometer.wave_length = Millers_Interferometer_Wave_Length;

	std::array<double,17> displacements;
	// turn the interferometer (by rotating velocity vector)
	for (int i=0;i<17;i++) {
		auto vi = rotate_z(v,rad(360.0/16*i)); // interferometer turns clockwise

		real t1 = phase_time(theory,vi,interferometer.light_path_1);
		real t2 = phase_time(theory,vi,interferometer.light_path_2);
		real dt = t1-t2;		
		displacements.at(i) = c*dt/interferometer.wave_length;
	}	

	if (invert) {
		for (auto& d : displacements)
			d = -d;
	}

	return displacements;
}






}//aether
