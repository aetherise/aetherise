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

	auto pL = dot_product(v,lp.p)/V;
	Vector3 x = v/V * pL;
	Vector3 y = lp.p-x;
	auto pT = length(y); // losing sign

	auto a = pL*V/sub_sqr(c/lp.n,V);
	auto b = (pL*pL+pT*pT)/sub_sqr(c/lp.n,V);
	return a + std::sqrt( sqr(a)+b );
}




real Aether::phase_delay(const Vector3& v,const Lightpath& lp) const
{
	const real c = Lightspeed;

	auto V = length(v); // velocity
	auto g = sqrt_1_minus_sqr(V/c);

	auto pL = dot_product(v,lp.p)/V;
	Vector3 x = v/V * pL;
	Vector3 y = lp.p-x;
	auto pT = length(y); // losing sign


	// anisotropic index of refraction
	auto Rd = (sqr(lp.n)-1)/(sqr(lp.n)+2);
	auto Rdg = Rd*g;
	auto nL = std::sqrt((1+2*Rdg)/(1-Rdg));

	auto s = length(lp.p);
	auto uL = pL/s*(c/nL);
	auto uT = pT/s*(c/lp.n);

	// velocity-addition formula		
	auto uL_ = (uL+V) / (1 + (uL/c)*(V/c));
	auto uT_ = (uT*g) / (1 + (uL/c)*(V/c));
	auto u_ = std::sqrt(sqr(uL_-V) + sqr(uT_));

	auto s_ = std::sqrt(sqr(pL*g) + sqr(pT));

	return s_/u_;
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



real phase_delay(const Theory& theory,const Vector3& v,const std::vector<Lightpath>& light_path)
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
	case Options::Theory::Aether:
		theory = make_unique<Aether>();
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

		real t1 = phase_delay(theory,vi,interferometer.light_path_1);
		real t2 = phase_delay(theory,vi,interferometer.light_path_2);
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
