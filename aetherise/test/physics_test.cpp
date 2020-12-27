
#include "test.h"

#include "../physics.h"
#include "../mathematics.h"

using namespace aether;


TEST(Celsius_Kelvis,test)
{
	{
		ASSERT_TRUE(in_Kelvin(0)==273.15);
		ASSERT_TRUE(in_Celsius(0)==-273.15);
		ASSERT_TRUE(in_Celsius(in_Kelvin(20)) == 20);
		ASSERT_TRUE(in_Kelvin(in_Celsius(20)) == 20);
	}
}





TEST(refractive_index_of_air,Table1)
{
	// Philip E. Ciddor: "Refractive index of air: new equations for the visible and near infrared", 1996
	// Table 1

	{
		auto n = refractive_index_of_air(80000,in_Kelvin(20),0.633,0,450);
		ASSERT_TRUE(approximate_delta((n-1)*1e8,21458,1));
	}

	{
		auto n = refractive_index_of_air(100000,in_Kelvin(20),0.633,0,450);
		ASSERT_TRUE(approximate_delta((n-1)*1e8,26824.4,1));
	}

	{
		auto n = refractive_index_of_air(120000,in_Kelvin(20),0.633,0,450);
		ASSERT_TRUE(approximate_delta((n-1)*1e8,32191.6,1));
	}

	{
		auto n = refractive_index_of_air(100000,in_Kelvin(10),0.633,0,450);
		ASSERT_TRUE(approximate_delta((n-1)*1e8,27774.7,1));
	}

	{
		auto n = refractive_index_of_air(100000,in_Kelvin(30),0.633,0,450);
		ASSERT_TRUE(approximate_delta((n-1)*1e8,25937.2,1));
	}
}



TEST(refractive_index_of_air,Table3)
{
	// Philip E. Ciddor: "Refractive index of air: new equations for the visible and near infrared", 1996
	// Table 3

	const auto prec=1;

	{
		auto n = refractive_index_of_air(80000,in_Kelvin(20),0.633,0.75,450);
		ASSERT_TRUE(approximate_delta((n-1)*1e8,21394,prec));
	}

	{
		auto n = refractive_index_of_air(120000,in_Kelvin(20),0.633,0.75,450);
		ASSERT_TRUE(approximate_delta((n-1)*1e8,32127.8,prec));
	}

	{
		auto n = refractive_index_of_air(80000,in_Kelvin(40),0.633,0.75,450);
		//ASSERT_TRUE(approximate_delta((n-1)*1e8,19996.5,prec)); // Ciddor value (probably typo)
		ASSERT_TRUE(approximate_delta((n-1)*1e8,19896.5,prec)); // calculator from NIST web page
	}

	{
		auto n = refractive_index_of_air(120000,in_Kelvin(40),0.633,0.75,450);
		ASSERT_TRUE(approximate_delta((n-1)*1e8,29941.8,prec));
	}

	{
		auto n = refractive_index_of_air(80000,in_Kelvin(50),0.633,1,450);
		ASSERT_TRUE(approximate_delta((n-1)*1e8,19058.4,prec));
	}

	{
		auto n = refractive_index_of_air(120000,in_Kelvin(50),0.633,1,450);
		ASSERT_TRUE(approximate_delta((n-1)*1e8,28792.4,prec));
	}
}



TEST(barometric_formula,test)
{
	{
		auto p = barometric_formula(0);
		ASSERT_TRUE(p == 101325);
	}

	{
		auto p = barometric_formula(1000);
		ASSERT_TRUE(approximate_delta(p,89876,1));
	}
}




int main(int ,char** )
{
	RUN_ALL_TESTS();

	return EXIT_SUCCESS;
}

