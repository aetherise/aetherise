

#include "test.h"

#include "../mathematics.h"
#include "../utils.h"
#include "../stdx.h"

#include <random>
#include <map>
#include <iomanip>

using namespace aether;



TEST(mean_value,test)
{
	{
		std::array<double,0> a;		
		auto m = mean_value(a.begin(),a.end());
		ASSERT_TRUE(std::isnan(m));		
	}

	{
		std::array<double,1> a {-1.1};
		auto m = mean_value(a.begin(),a.end());
		ASSERT_TRUE(m==-1.1);
	}

	{
		std::array<double,2> a {0,0};
		auto m = mean_value(a.begin(),a.end());
		ASSERT_TRUE(m==0);
	}

	{
		std::array<double,2> a {-1,1};
		auto m = mean_value(a.begin(),a.end());
		ASSERT_TRUE(m==0);
	}

	{
		std::array<double,3> a {0,3,6};
		auto m = mean_value(a.begin(),a.end());
		ASSERT_TRUE(m==3);
	}
}


TEST(mean_value,array)
{
	{
		std::array<double,3> a {1,3,2};
		auto m = mean_value(a);
		ASSERT_TRUE(m==2);
	}
}


TEST(median_value,vector)
{
	{
		std::vector<double> v;
		auto m = median_value(v);
		ASSERT_TRUE(std::isnan(m));
	}
	
	{
		std::vector<double> v {3};
		auto m = median_value(v);
		ASSERT_TRUE(m==3);
	}
	
	{
		std::vector<double> v {1,5};
		auto m = median_value(v);
		ASSERT_TRUE(m==3);
	}
	
	{
		std::vector<double> v {0,5,99};
		auto m = median_value(v);
		ASSERT_TRUE(m==5);
	}
	
	{
		std::vector<double> v {0,5,5,99};
		auto m = median_value(v);
		ASSERT_TRUE(m==5);
	}
	
	{
		std::vector<double> v {0,2,4,99};
		auto m = median_value(v);
		ASSERT_TRUE(m==3);
	}
	
	{
		std::vector<double> v {1,1,3,4,5,7,11};
		auto m = median_value(v);
		ASSERT_TRUE(m==4);
	}
}


TEST(rad,test)
{
	{
		auto r = rad(0);
		ASSERT_TRUE(r==0)
	}

	{
		auto r = rad(180);
		ASSERT_APPROX(r,AETHER_PI,0.000001);
	}

	{
		auto r = rad(360);
		ASSERT_APPROX(r,AETHER_2PI,0.000001);
	}
}



TEST(deg,test)
{
	{
		auto d = deg(0);
		ASSERT_TRUE(d==0)
	}

	{
		auto d = deg(AETHER_PI);
		ASSERT_APPROX(d,180.0,0.000001);
	}

	{
		auto d = deg(AETHER_2PI);
		ASSERT_APPROX(d,360.0,0.000001);
	}
}



TEST(rad_to_h,test)
{
	{
		auto h = rad_to_h(0);
		ASSERT_TRUE(h==0)
	}

	{
		auto h = rad_to_h(AETHER_PI);
		ASSERT_APPROX(h,12.0,0.000001);
	}

	{
		auto h = rad_to_h(AETHER_2PI);
		ASSERT_APPROX(h,24.0,0.000001);
	}
}



TEST(h_to_rad,test)
{
	{
		auto r = h_to_rad(0);
		ASSERT_TRUE(r==0)
	}

	{
		auto r = h_to_rad(12);
		ASSERT_APPROX(r,AETHER_PI,0.000001);
	}

	{
		auto r = h_to_rad(24);
		ASSERT_APPROX(r,AETHER_2PI,0.000001);
	}
}




TEST(KahanSum,test)
{
	{
		KahanSum<double> s;
		ASSERT_TRUE(s.sum()==0.0);
	}

	{
		KahanSum<double> s;
		s.add(1);
		s.add(1);
		ASSERT_TRUE(s.sum()==2.0);
	}

	{
		KahanSum<double> s;
		s.add(1);
		s.add(-1);
		ASSERT_TRUE(s.sum()==0.0);
	}

	{
		KahanSum<double> s;
		s.add(1e-17); // get lost
		s.add(1e+17);
		ASSERT_TRUE(s.sum()==1e17);
	}

	{
		KahanSum<double> s;
		s.add(1e-17); // still gets lost
		s.add(1e+17);
		s.add(-1e17);
		ASSERT_TRUE(s.sum()==0.0);
	}
	
	{
		KahanSum<double> s;
		s.add(1);
		s.add(1e-16);								
		ASSERT_TRUE(s.sum()==1.0);
	}
	
	{
		KahanSum<double> s;
		s.add(1);
		s.add(1e-16);			
		s.add(1e-16);			
		s.add(1e-16);			
		
		ASSERT_TRUE(s.sum()>1.0);
	}

	{
		KahanSum<double> s;
		s+=1;
		s+=1;
		ASSERT_TRUE(s.sum()==2.0);
	}

	
}




TEST(Vector3,test)
{
	{
		Vector3 v {1,2,3};
		ASSERT_TRUE(v.x==1 && v.y==2 && v.z==3);
	}
}

TEST(Vector3,operators)
{
	{
		Vector3 a {1,2,3};
		Vector3 b {-1,-1,-1};

		Vector3 c = a+b;
		ASSERT_TRUE(c.x==0 && c.y==1 && c.z==2);
	}


	{
		Vector3 a {1,2,3};
		Vector3 b {-1,-1,-1};

		Vector3 c = a-b;
		ASSERT_TRUE(c.x==2 && c.y==3 && c.z==4);
	}

	{
		Vector3 a {1,2,3};

		Vector3 c = a*3;
		ASSERT_TRUE(c.x==3 && c.y==6 && c.z==9);
	}

	{
		Vector3 a {1,2,3};

		Vector3 c = a/2;
		ASSERT_TRUE(c.x==0.5 && c.y==1 && c.z==1.5);
	}

}


TEST(dot_product,test)
{
	{
		Vector3 a {1,1,0};
		Vector3 b {-1,1,0};

		real p = dot_product(a,b);
		ASSERT_TRUE(p==0);
	}

	{
		Vector3 a {0,1,1};
		Vector3 b {0,-1,1};

		real p = dot_product(a,b);
		ASSERT_TRUE(p==0);
	}


	{
		Vector3 a {1,2,3};
		Vector3 b {-1,-1,-1};

		real p = dot_product(a,b);
		ASSERT_TRUE(p==-6);
	}
}


TEST(length,test)
{
	{
		Vector3 v {0,0,0};
		auto l = length(v);
		ASSERT_TRUE(l==0);
	}

	{
		Vector3 v {1,0,0};
		auto l = length(v);
		ASSERT_TRUE(l==1);
	}

	{
		Vector3 v {0,1,0};
		auto l = length(v);
		ASSERT_TRUE(l==1);
	}

	{
		Vector3 v {0,0,1};
		auto l = length(v);
		ASSERT_TRUE(l==1);
	}

	{
		Vector3 v {1,0,1};
		auto l = length(v);
		ASSERT_APPROX(l,std::sqrt(2),0.000001);
	}

	{
		Vector3 v {0,3,4};
		auto l = length(v);
		ASSERT_TRUE(l==5);
	}

	{
		Vector3 v {2,2,2};
		auto l = length(v);
		ASSERT_APPROX(l,2*std::sqrt(3),0.000001);
	}
}


TEST(rotate_x,test)
{
	{
		Vector3 v {1,0,0};
		Vector3 w = rotate_x(v,rad(90));
		ASSERT_TRUE(approximate(w.x,1,0.000001));
		ASSERT_TRUE(approximate(w.y,0,0.000001));
		ASSERT_TRUE(approximate(w.z,0,0.000001));
	}

	{
		Vector3 v {0,1,0};
		Vector3 w = rotate_x(v,rad(90));
		ASSERT_TRUE(approximate(w.x,0,0.000001));
		ASSERT_TRUE(approximate(w.y,0,0.000001));
		ASSERT_TRUE(approximate(w.z,1,0.000001));
	}


}



TEST(rotate_y,test)
{
	{
		Vector3 v {0,1,0};
		Vector3 w = rotate_y(v,rad(90));
		ASSERT_TRUE(approximate(w.x,0,0.000001));
		ASSERT_TRUE(approximate(w.y,1,0.000001));
		ASSERT_TRUE(approximate(w.z,0,0.000001));
	}

	{
		Vector3 v {0,0,1};
		Vector3 w = rotate_y(v,rad(90));
		ASSERT_TRUE(approximate(w.x,-1,0.000001)); // correct?
		ASSERT_TRUE(approximate(w.y,0,0.000001));
		ASSERT_TRUE(approximate(w.z,0,0.000001));
	}


}


TEST(rotate_z,test)
{
	{
		Vector3 v {0,0,1};
		Vector3 w = rotate_z(v,rad(90));
		ASSERT_TRUE(approximate(w.x,0,0.000001));
		ASSERT_TRUE(approximate(w.y,0,0.000001));
		ASSERT_TRUE(approximate(w.z,1,0.000001));
	}

	{
		Vector3 v {1,0,0};
		Vector3 w = rotate_z(v,rad(90));
		ASSERT_TRUE(approximate(w.x,0,0.000001));
		ASSERT_TRUE(approximate(w.y,1,0.000001));
		ASSERT_TRUE(approximate(w.z,0,0.000001));
	}
}


TEST(sqr,test)
{
	{
		ASSERT_TRUE(sqr(1)==1);
	}

	{
		ASSERT_TRUE(sqr(-1)==1);
	}


	{
		ASSERT_TRUE(sqr(2)==4);
	}


	{
		ASSERT_APPROX(sqr(std::sqrt(3)),3,0.000001);
	}

}



TEST(sub_sqr,test)
{
	{
		auto y = sub_sqr(0,0);
		ASSERT_TRUE(y==0);
	}

	{
		auto y = sub_sqr(1,1);
		ASSERT_TRUE(y==0);
	}

	{
		auto y = sub_sqr(1,2);
		ASSERT_TRUE(y==-3);
	}

	{
		auto y = sub_sqr(2,1);
		ASSERT_TRUE(y==3);
	}

}



TEST(_1_minus_sqr,test)
{
	{
		auto y = _1_minus_sqr(0);
		ASSERT_TRUE(y==1);
	}

	{
		auto y = _1_minus_sqr(1);
		ASSERT_TRUE(y==0);
	}

	{
		auto y = _1_minus_sqr(-1);
		ASSERT_TRUE(y==0);
	}


	{
		auto y = _1_minus_sqr(2);
		ASSERT_TRUE(y==-3);
	}

	{
		auto y = _1_minus_sqr(std::sqrt(2));
		ASSERT_APPROX(y,-1,0.000001);
	}

}


TEST(sqrt_1_minus_sqr,test)
{
	{
		auto y = sqrt_1_minus_sqr(0);
		ASSERT_TRUE(y==1);
	}

	{
		auto y = sqrt_1_minus_sqr(1);
		ASSERT_TRUE(y==0);
	}

	{
		auto y = sqrt_1_minus_sqr(-1);
		ASSERT_TRUE(y==0);
	}


	{
		auto y = sqrt_1_minus_sqr(1/std::sqrt(2));
		ASSERT_APPROX(y,std::sqrt(0.5),0.000001);
	}
}


TEST(period_360,test)
{
	{
		auto p = period_360(0);
		ASSERT_TRUE(p==0);
	}

	{
		auto p = period_360(-180);
		ASSERT_TRUE(p==180);
	}

	{
		auto p = period_360(180);
		ASSERT_TRUE(p==180);
	}

	{
		auto p = period_360(360);
		ASSERT_TRUE(p==0);
	}

	{
		auto p = period_360(-360);
		ASSERT_TRUE(p==0);
	}

	{
		auto p = period_360(360+180);
		ASSERT_TRUE(p==180);
	}

	{
		auto p = period_360(2*360+180);
		ASSERT_TRUE(p==180);
	}

	{
		auto p = period_360(-360-180);
		ASSERT_TRUE(p==180);
	}

	{
		auto p = period_360(-2*360-180);
		ASSERT_TRUE(p==180);
	}

}




TEST(period_2pi,test)
{
	{
		auto p = period_2pi(0);
		ASSERT_TRUE(p==0);
	}

	{
		auto p = period_2pi(-AETHER_PI);
		ASSERT_APPROX(p,AETHER_PI,0.000001);
	}

	{
		auto p = period_2pi(AETHER_PI);
		ASSERT_APPROX(p,AETHER_PI,0.000001);
	}

	{
		auto p = period_2pi(AETHER_2PI);
		ASSERT_APPROX(p,0,0.000001);
	}

	{
		auto p = period_2pi(AETHER_2PI);
		ASSERT_APPROX(p,0,0.000001);
	}

	{
		auto p = period_2pi(AETHER_2PI+AETHER_PI);
		ASSERT_APPROX(p,AETHER_PI,0.000001);
	}

	{
		auto p = period_2pi(2*AETHER_2PI+AETHER_PI);
		ASSERT_APPROX(p,AETHER_PI,0.000001);
	}

	{
		auto p = period_2pi(-AETHER_2PI-AETHER_PI);
		ASSERT_APPROX(p,AETHER_PI,0.000001);
	}

	{
		auto p = period_2pi(-2*AETHER_2PI-AETHER_PI);
		ASSERT_APPROX(p,AETHER_PI,0.000001);
	}

}



TEST(period_24,test)
{
	{
		auto p = period_24(0);
		ASSERT_TRUE(p==0);
	}

	{
		auto p = period_24(-12);
		ASSERT_TRUE(p==12);
	}

	{
		auto p = period_24(12);
		ASSERT_TRUE(p==12);
	}

	{
		auto p = period_24(24);
		ASSERT_TRUE(p==0);
	}

	{
		auto p = period_24(-24);
		ASSERT_TRUE(p==0);
	}

	{
		auto p = period_24(24+12);
		ASSERT_TRUE(p==12);
	}

	{
		auto p = period_24(2*24+12);
		ASSERT_TRUE(p==12);
	}

	{
		auto p = period_24(-24-12);
		ASSERT_TRUE(p==12);
	}

	{
		auto p = period_24(-2*24-12);
		ASSERT_TRUE(p==12);
	}

}


TEST(polar3,test)
{
	{
		Vector3 v {0,0,0};
		auto p = polar3(v);
		ASSERT_TRUE(p.r==0);
	}

	{
		Vector3 v {1,0,0};
		auto p = polar3(v);

		ASSERT_APPROX(p.r,1,0.000001);
		ASSERT_APPROX(p.th,rad(90),0.000001);
		ASSERT_APPROX(p.ph,0,0.000001);
	}
	
	{
		Vector3 v {-1,0,0};
		auto p = polar3(v);

		ASSERT_APPROX(p.r,1,0.000001);
		ASSERT_APPROX(p.th,rad(90),0.000001);
		ASSERT_APPROX(p.ph,rad(180),0.000001);
	}
	

	{
		Vector3 v {0,1,0};
		auto p = polar3(v);

		ASSERT_APPROX(p.r,1,0.000001);
		ASSERT_APPROX(p.th,rad(90),0.000001);
		ASSERT_APPROX(p.ph,rad(90),0.000001);
	}
	
	{
		Vector3 v {0,-1,0};
		auto p = polar3(v);

		ASSERT_APPROX(p.r,1,0.000001);
		ASSERT_APPROX(p.th,rad(90),0.000001);
		ASSERT_APPROX(period_2pi(p.ph),rad(270),0.000001);
	}

	{
		Vector3 v {0,0,2};
		auto p = polar3(v);

		ASSERT_APPROX(p.r,2,0.000001);
		ASSERT_APPROX(p.th,rad(0),0.000001);
		//ASSERT_APPROX(p.ph,rad(0),0.000001); // undefined
	}

	{
		Vector3 v {0,0,-2};
		auto p = polar3(v);

		ASSERT_APPROX(p.r,2,0.000001);
		ASSERT_APPROX(p.th,rad(180),0.000001);
		//ASSERT_APPROX(p.ph,rad(0),0.000001); // undefined
	}
}



TEST(vector3,test)
{
	{
		for (int phi=0;phi<=360;phi+=15) {
			for (int theta=15;theta<=165;theta+=15) {
				auto v = vector3({2,rad(theta),rad(phi)});
				auto p = polar3(v);
				ASSERT_APPROX(p.r,2,1e-6);
				ASSERT_APPROX(p.th,rad(theta),1e-6);
				ASSERT_TRUE(approximate_delta_periodic(p.ph,rad(phi),1e-5,AETHER_2PI));
			}
		}
	}
}




TEST(integrate,test)
{
	{

		auto s = integrate(0,0,[](real ){
			return 0;
		},3);

		ASSERT_TRUE(s==0);
	}

	{
		auto s = integrate(1,2,[](real ){
			return 0;
		},3);

		ASSERT_TRUE(s==0);
	}

	{
		auto s = integrate(1,2,[](real ){
			return 1;
		},3);

		ASSERT_APPROX(s,1,0.001);
	}

	{
		auto s = integrate(1,1,[](real ){
			return 1;
		},3);

		ASSERT_APPROX(s,0,0.001);
	}


	{
		auto s = integrate(2,1,[](real ){
			return 1;
		},3);

		ASSERT_APPROX(s,-1,0.001);
	}
}



TEST(integrate,infnan)
{
	{
		auto s = integrate(-1,1,[](real ){
			return INFINITY;
		},3);

		ASSERT_TRUE(std::isnan(s));
	}

	{
		auto s = integrate(-1,1,[](real ){
			return NAN;
		},3);

		ASSERT_TRUE(std::isnan(s));
	}
}




TEST(chi_squared_test,test)
{
	{
		std::array<double,16> data {};
		std::array<double,16> uncertainty {};
		std::array<double,16> theory {};
		auto chi = chi_squared_test(data,uncertainty,theory);
		ASSERT_TRUE(std::isnan(chi));
	}

	{
		std::array<double,16> data {};
		std::array<double,16> uncertainty {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
		std::array<double,16> theory {};
		auto chi = chi_squared_test(data,uncertainty,theory);
		ASSERT_TRUE(chi==0);
	}

	{
		std::array<double,16> data {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
		std::array<double,16> uncertainty {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
		std::array<double,16> theory {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
		auto chi = chi_squared_test(data,uncertainty,theory);
		ASSERT_TRUE(chi==0);
	}

	{
		std::array<double,16> data        {4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4};
		std::array<double,16> uncertainty {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
		std::array<double,16> theory      {2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2};
		auto chi = chi_squared_test(data,uncertainty,theory);
		ASSERT_TRUE(chi==16*4);
	}

}



TEST(integrate,examples)
{

	{
		// Example from Papula 4. Edition
		auto s = integrate(0,1,[](real x){
			return std::exp(-x*x);
		},3);

		ASSERT_APPROX(s,0.746824,0.000001);
	}

	

	{
		auto s = integrate(-1,1,[](real x){
			return std::sqrt(1-x*x); // half circle
		},8);
		ASSERT_APPROX(s,AETHER_PI_2,0.1);
	}

	
	{
		// https://de.wikipedia.org/wiki/Romberg-Integration#cite_note-1
		int n=8;
		auto s = integrate(0,AETHER_2PI,[&](real x){
			return std::cos(std::pow(2,n)*x);
		},8);

		// not 0
		ASSERT_FALSE(test::approximate(s,0,0.0001));
	}

	{
		// https://de.wikipedia.org/wiki/Romberg-Integration#cite_note-1
		int n=8;
		auto s = integrate(0,AETHER_2PI,[&](real x){
			return std::cos(std::pow(2,n)*x);
		},n+4);

		// gets 0
		ASSERT_APPROX(s,0,0.0001);
	}

}



TEST(chi_square_cdf,test)
{
	{
		double chi = 9;
		int n = 14; // even

		auto s1 = integrate(0.,chi,[n](real x){
			return chi_square_pdf(x,n);
		},10);

		auto s2 = chi_square_cdf(chi,n);

		ASSERT_APPROX(s1,s2,0.00001);
	}

	{
		double chi = 9;
		int n = 15; // odd

		auto s1 = integrate(0.,chi,[n](real x){
			return chi_square_pdf(x,n);
		},10);

		auto s2 = chi_square_cdf(chi,n);

		ASSERT_APPROX(s1,s2,0.00001);
	}

	{
		auto p = 1-chi_square_cdf(9.34,10);
		ASSERT_APPROX(p,0.5,0.01);
	}

	{
		auto p = 1-chi_square_cdf(18.31,10);
		ASSERT_APPROX(p,0.05,0.01);
	}
	
	// high n
	{
		auto p = 1-chi_square_cdf(500,500);
		ASSERT_APPROX(p,0.5,0.1);	
	}
	{
		auto p = 1-chi_square_cdf(1000,1000);
		ASSERT_APPROX(p,0.5,0.1);	
	}
	{
		auto p = 1-chi_square_cdf(1500,1500);
		ASSERT_APPROX(p,0.5,0.1);	
	}
}



TEST(norm_cdf,test)
{
	{
		double m = 0;
		double s = 1;
		auto s1 = integrate(-10,10,[m,s](real x){
			return normal_pdf(m,s,x);
		},10);

		auto s2 = normal_cdf(m,s,10);

		ASSERT_APPROX(s1,s2,0.00001);
		ASSERT_APPROX(s1,1.,0.00001);
	}
	
	{
		double m = 2;
		double s = 1;
		auto s1 = integrate(-8,12,[m,s](real x){
			return normal_pdf(m,s,x);
		},10);

		auto s2 = normal_cdf(m,s,10);

		ASSERT_APPROX(s1,s2,0.00001);
		ASSERT_APPROX(s1,1.,0.00001);
	}
	
	{
		double m = 2;
		double s = 1;
		auto s1 = integrate(1,3,[m,s](real x){
			return normal_pdf(m,s,x);
		},10);

		auto s2 = normal_cdf(m,s,3)-normal_cdf(m,s,1);

		ASSERT_APPROX(s1,s2,0.00001);
		ASSERT_APPROX(s1,0.6827,0.001);
	}
}



TEST(test_for_normality,test)
{
	// Stephens
	//const double q_10 = ADTestStephensQuantiles.q_10;
	//const double q_5 = ADTestStephensQuantiles.q_5;

	// D'Agostino
	//const double q_10 = 0.249; // 75%
	const double q_10 = ADTestDAgostinoQuantiles.q_10;
	const double q_5 = ADTestDAgostinoQuantiles.q_5;

	double percentnn = 0.02;
	double syserr = 0.2;

	{
		// M. A. Stephens (1974): EDF Statistics for Goodness of Fit and Some Comparisons
		// Journal of the American Statistical Association, 69:347, 730-737
		auto A = test_for_normality({148,154,158,160,161,162,166,170,182,195,236},ADTestType::Stephens);
		ASSERT_APPROX(A,1.095,0.001);
	}

	{
		std::vector<double> x;
		std::mt19937 engine(create_random_seed());
		std::normal_distribution<double> dist(0,0.02);

		const int n = 10000;
		int nnd = 0;
		for (int i=0;i<n;i++) {
			x.clear();
			for (int k=0;k<20;k++) {
				x.push_back(dist(engine));
			}
			std::sort(x.begin(),x.end());
			auto A = test_for_normality(x);
			if (A > q_5)
				nnd++;
		}

		std::cout << "not normal distributed (5%): " << 100.0*nnd/n << "%\n";
		//ASSERT_APPROX(nnd,500,0.5);
	}

	{
		std::vector<double> x;
		std::mt19937 engine(create_random_seed());
		std::normal_distribution<double> dist(0,0.02);

		const int n = 10000;
		int nnd = 0;
		for (int i=0;i<n;i++) {
			x.clear();
			for (int k=0;k<20;k++) {
				x.push_back(dist(engine));
			}
			std::sort(x.begin(),x.end());
			auto A = test_for_normality(x);
			if (A > q_10)
				nnd++;
		}

		std::cout << "not normal distributed (10%): " << 100.0*nnd/n << "%\n";
		//ASSERT_APPROX(nnd,1000,0.5);
	}

	{
		std::vector<double> x;
		std::mt19937 engine(create_random_seed());
		std::normal_distribution<double> dist(0,0.02);
		std::uniform_real_distribution<double> dist_nnd(-0.04,0.04);

		const int n = 10000;
		int nnd = 0;
		for (int i=0;i<n;i++) {
			x.clear();
			if (i<n-n*percentnn) { // 2% not normal
				// normal
				for (int k=0;k<20;k++) {
					x.push_back(dist(engine));
				}
			}
			else {
				if (false) {
					// uniform
					for (int k=0;k<20;k++) {
						// will be fully rejected at 75% level
						// will not be detected at 5% level
						x.push_back(dist_nnd(engine));
					}
				}
				else {
					// double normal: systematic error change
					for (int k=0;k<10;k++) {
						x.push_back(dist(engine));
					}
					for (int k=0;k<10;k++) {
						x.push_back(dist(engine)+syserr);
					}
				}
			}
			std::sort(x.begin(),x.end());
			auto A = test_for_normality(x);
			if (A > q_5)
				nnd++;
		}

		std::cout << "not normal distributed (7%): " << 100.0*nnd/n << "%\n";
		//ASSERT_APPROX(nnd,700,0.2);
	}

	{
		std::vector<double> x;
		std::mt19937 engine(create_random_seed());
		std::normal_distribution<double> dist(0,0.02);
		std::uniform_real_distribution<double> dist_nnd(-0.04,0.04);

		const int n = 10000;
		int nnd = 0;
		for (int i=0;i<n;i++) {
			x.clear();
			if (i<n-n*percentnn) { // 2% not normal
				// normal
				for (int k=0;k<20;k++) {
					x.push_back(dist(engine));
				}
			}
			else {
				if (false) {
					// uniform
					for (int k=0;k<20;k++) {
						// will be fully rejected at 75% level
						// will not be detected at 5% level
						x.push_back(dist_nnd(engine));
					}
				}
				else {
					// double normal: systematic error change
					for (int k=0;k<10;k++) {
						x.push_back(dist(engine));
					}
					for (int k=0;k<10;k++) {
						x.push_back(dist(engine)+syserr);
					}
				}
			}
			std::sort(x.begin(),x.end());
			auto A = test_for_normality(x);
			if (A > q_10) // 10%
				nnd++;
		}

		std::cout << "not normal distributed (12%): " << 100.0*nnd/n << "%\n";
		//ASSERT_APPROX(nnd,1200,0.4);
	}
}



TEST(confidence_interval,test)
{
	{
		int k = 50;
		int n = 100;
		ConfidenceInterval ci = confidence_interval(k,n);
		
		auto p_ = (k+2.)/(n+4.); // approximation
		auto u_ = 2*std::sqrt((p_-p_*p_)/(n+4)); // approximation
		ASSERT_APPROX(ci.p,p_,0.01);
		ASSERT_APPROX(ci.u,u_,0.02);
		
	}
	
	{
		int k = 25;
		int n = 100;
		ConfidenceInterval ci1 = confidence_interval(k,n);
		ConfidenceInterval ci2 = confidence_interval(n-k,n);
				
		ASSERT_APPROX(ci1.p,1-ci2.p,0.000001);		
		ASSERT_APPROX(ci1.u,ci2.u,0.000001);		
		
	}
}



TEST(p_value_A,test)
{
	const auto q = test_quantiles(ADTestType::DAgostino);

	{
		auto p = p_value_A(0.188);
		ASSERT_APPROX(p,0.90,0.1);
	}

	{
		auto p = p_value_A(0.249);
		ASSERT_APPROX(p,0.75,0.1);
	}

	{
		auto p = p_value_A(q.q_50);
		ASSERT_APPROX(p,0.5,0.1);
	}

	{
		auto p = p_value_A(q.q_25);
		ASSERT_APPROX(p,0.25,0.1);
	}

	{
		auto p = p_value_A(q.q_15);
		ASSERT_APPROX(p,0.15,0.1);
	}

	{
		auto p = p_value_A(q.q_10);
		ASSERT_APPROX(p,0.10,0.1);
	}

	{
		auto p = p_value_A(q.q_5);
		ASSERT_APPROX(p,0.05,0.1);
	}

	{
		auto p = p_value_A(q.q_1);
		ASSERT_APPROX(p,0.01,0.1);
	}
}



TEST(sgn,test)
{
	{
		ASSERT_TRUE(sgn(0)==0);
		ASSERT_TRUE(sgn(-1)==-1);
		ASSERT_TRUE(sgn(1)==1);
		ASSERT_TRUE(sgn(-13)==-1);
		ASSERT_TRUE(sgn(7)==1);		
	}
	
	{
		ASSERT_TRUE(sgn(0.0)==0);
		ASSERT_TRUE(sgn(-1.0)==-1);
		ASSERT_TRUE(sgn(1.0)==1);
		ASSERT_TRUE(sgn(-13.1)==-1);
		ASSERT_TRUE(sgn(7.1)==1);		
	}
}



TEST(sample_variance,test)
{
	{
		std::vector<double> data;		
		auto s = sample_variance(data.begin(),data.end(),1.);
		ASSERT_TRUE(std::isnan(s));
	}
	
	{
		std::vector<double> data {1,1,1,1,1};
		auto mean = mean_value(data.begin(),data.end());
		auto s = sample_variance(data.begin(),data.end(),mean);
		ASSERT_APPROX(s,0,0.00001);
	}
	
	{
		std::vector<double> data {3,1,1,3};
		auto mean = mean_value(data.begin(),data.end());
		auto s = sample_variance(data.begin(),data.end(),mean);
		ASSERT_APPROX(s,4.0/3.0,0.00001);
	}
}



TEST(sample_variance,array)
{
	{
		std::array<double,4> data {3,1,1,3};
		auto mean = mean_value(data);
		auto s = sample_variance(data,mean);
		ASSERT_APPROX(s,4.0/3.0,0.00001);
	}
}



TEST(sample_standard_deviation,test)
{
	{
		std::vector<double> data;		
		auto s = sample_standard_deviation(data.begin(),data.end(),1.);
		ASSERT_TRUE(std::isnan(s));
	}
	
	{
		std::vector<double> data {1,1,1,1,1};
		auto mean = mean_value(data.begin(),data.end());
		auto s = sample_standard_deviation(data.begin(),data.end(),mean);
		ASSERT_APPROX(s,0,0.00001);
	}
	
	{
		std::vector<double> data {3,1,1,3};
		auto mean = mean_value(data.begin(),data.end());
		auto s = sample_standard_deviation(data.begin(),data.end(),mean);
		ASSERT_APPROX(s,std::sqrt(4.0/3.0),0.00001);
	}
}



TEST(sample_standard_deviation,array)
{	
	{
		std::array<double,4> data {3,1,1,3};
		auto mean = mean_value(data);
		auto s = sample_standard_deviation(data,mean);
		ASSERT_APPROX(s,std::sqrt(4.0/3.0),0.00001);
	}
}





TEST(bisection,test)
{
	{
		double x=-2;
		bisection(-1.0,1.0,[&x](double ,double p,double ){
			x=p;
			return 0;
		});
		ASSERT_APPROX(x,0,0.000000001);
	}
	
	{
		double x=-1;
		bisection(0.0,0.0,[&x](double ,double p,double ){
			x=p;
			return 1.0;
		});
		ASSERT_TRUE(x==-1); // no call
	}
	
	{
		double x=-2;
		bisection(-1.0,1.0,[&x](double ,double p,double ){
			x=p;
			return 1.0;
		});
		ASSERT_APPROX(x,1,0.000000001);
	}
	
	{
		double x=-2;
		bisection(1.0,-1.0,[&x](double ,double p,double ){
			x=p;
			return 1.0;
		});
		ASSERT_APPROX(x,-1,0.000000001);
	}
	
	
	{
		double x=-1;
		bisection(0.0,0.0,[&x](double ,double p,double ){
			x=p;
			return -1.0;
		});
		ASSERT_TRUE(x==-1); // no call
	}
	
	{
		double x=-2;
		bisection(-1.0,1.0,[&x](double ,double p,double ){
			x=p;
			return -1.0;
		});
		ASSERT_APPROX(x,-1,0.000000001);
	}
	
	{
		double x=-2;
		bisection(1.0,-1.0,[&x](double ,double p,double ){
			x=p;
			return -1.0;
		});
		ASSERT_APPROX(x,1,0.000000001);
	}
}




TEST(derive,polynomials)
{
	{
		double h = 0.25;
		for (double x=-2;x<=2;x+=0.1) {

			double y_x = derive(x,[](double x){
				return 2*x*x-3*x+5;
			},h);

			auto y_ = [](double x){
				return 4*x-3;
			};

			ASSERT_APPROX(y_x,y_(x),0.000001);
		}
	}
}



TEST(derive,sin)
{
	{
		double x = 0;
		double h = 0.000001;
		auto y_ = derive(x,[](double x){
			return sin(x);
		},h);
		ASSERT_APPROX(y_,cos(x),1e-6);
	}

	{
		double x = 1;
		double h = 0.000001;
		auto y_ = derive(x,[](double x){
			return sin(x);
		},h);
		ASSERT_APPROX(y_,cos(x),1e-6);
	}

	{
		double x = 2;
		double h = 0.000001;
		auto y_ = derive(x,[](double x){
			return sin(x);
		},h);
		ASSERT_APPROX(y_,cos(x),1e-6);
	}

	{
		double x = 3;
		double h = 0.000001;
		auto y_ = derive(x,[](double x){
			return sin(x);
		},h);
		ASSERT_APPROX(y_,cos(x),1e-6);
	}

	{
		double x = 4;
		double h = 0.000001;
		auto y_ = derive(x,[](double x){
			return sin(x);
		},h);
		ASSERT_APPROX(y_,cos(x),1e-6);
	}
	
}



TEST(grad,test)
{
	{
		double h = 0.001;
		std::array<double,2> args {0,0};
		
		auto f = [](const std::array<double,2>& x) {
			return x[0]*x[1];
		};
		auto g = grad(f,args,h);
		ASSERT_APPROX(g[0],0,1e-9);
		ASSERT_APPROX(g[1],0,1e-9);
	}
	
	{
		double h = 0.001;
		std::array<double,2> args {1,0};
		
		auto f = [](const std::array<double,2>& x) {
			return x[0]*x[1];
		};
		auto g = grad(f,args,h);
		ASSERT_APPROX(g[0],0,1e-9);
		ASSERT_APPROX(g[1],1,1e-9);
	}
	
	{
		double h = 0.001;
		std::array<double,2> args {0,-2};
		
		auto f = [](const std::array<double,2>& x) {
			return x[0]*x[1];
		};
		auto g = grad(f,args,h);
		ASSERT_APPROX(g[0],-2,1e-9);
		ASSERT_APPROX(g[1],0,1e-9);
	}
}



TEST(periodic_distance,test)
{
	{
		ASSERT_APPROX(periodic_distance(0,2,24), 2, 0.000001);
		ASSERT_APPROX(periodic_distance(2,0,24), 2, 0.000001);
		ASSERT_APPROX(periodic_distance(23,1,24), 2, 0.000001);
		ASSERT_APPROX(periodic_distance(1,23,24), 2, 0.000001);
		ASSERT_APPROX(periodic_distance(22,24,24), 2, 0.000001);
		ASSERT_APPROX(periodic_distance(24,22,24), 2, 0.000001);
		ASSERT_APPROX(periodic_distance(0,12,24), 12, 0.000001);
		ASSERT_APPROX(periodic_distance(12,0,24), 12, 0.000001);
		ASSERT_APPROX(periodic_distance(12,24,24), 12, 0.000001);
		ASSERT_APPROX(periodic_distance(24,12,24), 12, 0.000001);
		ASSERT_APPROX(periodic_distance(0,24,24), 0, 0.000001);
	}
}





TEST(DFT_analyze,test)
{
	const int samples = 36;
	std::vector<double> a;
	for (int i=0;i<samples;i++) {
		a.push_back(2*std::cos(AETHER_2PI/samples*(i)) + std::cos(2*AETHER_2PI/samples*(i)+2) + std::cos(3*AETHER_2PI/samples*(i)+5));
	}
	
	{		
		const int k = 1;
		auto y = DFT_analyze(k,a.begin(),a.end());		
		auto amplitude = std::abs(y);
		auto phase =std::arg(y);		
		std::cout << "amplitude = " << amplitude << "\n";
		std::cout << "phase = " << phase << "\n";
		ASSERT_APPROX(amplitude,2,0.001);		
		ASSERT_APPROX(period_2pi(phase),0,0.001);		
	}
	
	{		
		const int k = 2;
		auto y = DFT_analyze(k,a.begin(),a.end());		
		auto amplitude = std::abs(y);
		auto phase = std::arg(y);		
		std::cout << "amplitude = " << amplitude << "\n";
		std::cout << "phase = " << phase << "\n";
		ASSERT_APPROX(amplitude,1,0.001);		
		ASSERT_APPROX(period_2pi(phase),2,0.001);		
	}
	
	{		
		const int k = 3;
		auto y = DFT_analyze(k,a.begin(),a.end());		
		auto amplitude = std::abs(y);
		auto phase = std::arg(y);		
		std::cout << "amplitude = " << amplitude << "\n";
		std::cout << "phase = " << phase << "\n";
		ASSERT_APPROX(amplitude,1,0.001);
		ASSERT_APPROX(period_2pi(phase),5,0.001);				
	}
	
	{		
		const int k = 4;
		auto y = DFT_analyze(k,a.begin(),a.end());		
		auto amplitude = std::abs(y);
		auto phase = std::arg(y);		
		std::cout << "amplitude = " << amplitude << "\n";
		std::cout << "phase = " << phase << "\n";
		ASSERT_APPROX(amplitude,0,0.001);		
	}	
}





TEST(DFTGoertzel,test)
{
	const int sample_rate = 36;
	const int samples = sample_rate*2;
	std::vector<double> a;
	for (int i=0;i<samples;i++) {
		a.push_back(2*std::cos(AETHER_2PI/sample_rate*(i)) + std::cos(2*AETHER_2PI/sample_rate*(i)+2) + std::cos(3*AETHER_2PI/sample_rate*(i)+5));
	}
	
	{		
		const double f = 1;		
		DFTGoertzel dft({f},sample_rate);
		dft.analyze(a.begin(),a.end());		
		auto y = dft.result().at(f);
		
		auto amplitude = std::abs(y);
		auto phase =std::arg(y);		
		std::cout << "amplitude = " << amplitude << "\n";
		std::cout << "phase = " << phase << "\n";
		ASSERT_APPROX(amplitude,2,0.001);		
		ASSERT_APPROX(period_2pi(phase),0,0.001);		
	}
	
	{		
		const double f = 2;
		DFTGoertzel dft({f},sample_rate);
		dft.analyze(a.begin(),a.end());		
		auto y = dft.result().at(f);
		
		auto amplitude = std::abs(y);
		auto phase = std::arg(y);		
		std::cout << "amplitude = " << amplitude << "\n";
		std::cout << "phase = " << phase << "\n";
		ASSERT_APPROX(amplitude,1,0.001);		
		ASSERT_APPROX(period_2pi(phase),2,0.001);		
	}
	
	{		
		const double f = 3;
		DFTGoertzel dft({f},sample_rate);
		dft.analyze(a.begin(),a.end());		
		auto y = dft.result().at(f);
		
		auto amplitude = std::abs(y);
		auto phase = std::arg(y);		
		std::cout << "amplitude = " << amplitude << "\n";
		std::cout << "phase = " << phase << "\n";
		ASSERT_APPROX(amplitude,1,0.001);
		ASSERT_APPROX(period_2pi(phase),5,0.001);				
	}
	
	{		
		const double f = 4;
		DFTGoertzel dft({f},sample_rate);
		dft.analyze(a.begin(),a.end());		
		auto y = dft.result().at(f);
		
		auto amplitude = std::abs(y);
		auto phase = std::arg(y);		
		std::cout << "amplitude = " << amplitude << "\n";
		std::cout << "phase = " << phase << "\n";
		ASSERT_APPROX(amplitude,0,0.001);		
	}	
}




int main(int ,char** )
{
	RUN_ALL_TESTS();

	return EXIT_SUCCESS;
}
