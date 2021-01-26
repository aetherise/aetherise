
#include "test.h"

#include "../utils.h"

using namespace aether;


TEST(split,with_quote)
{
	{
		std::string s = "\";\"";
		auto v = split(s,';');
		ASSERT_TRUE(v.size()==1);
		ASSERT_TRUE(v.at(0)=="\";\"");
	}

	{
		std::string s = ";\";\"";
		auto v = split(s,';');
		ASSERT_TRUE(v.size()==2);
		ASSERT_TRUE(v.at(1)=="\";\"");
	}

	{
		std::string s = ";\";\";";
		auto v = split(s,';');
		ASSERT_TRUE(v.size()==3);
		ASSERT_TRUE(v.at(1)=="\";\"");
	}
}



TEST(split,with_brackets)
{
	{
		std::string s = "";
		auto tokens = split(s,';','(',')');
		ASSERT_TRUE(tokens.size()==1u);
		ASSERT_TRUE(tokens.at(0).empty());
	}
	
	{
		std::string s = ",";
		auto tokens = split(s,',','(',')');
		ASSERT_TRUE(tokens.size()==2u);
		ASSERT_TRUE(tokens.at(0).empty());
		ASSERT_TRUE(tokens.at(1).empty());
	}
	
	{
		std::string s = "a,";
		auto tokens = split(s,',','(',')');
		ASSERT_TRUE(tokens.size()==2u);
		ASSERT_TRUE(tokens.at(0) == "a");
		ASSERT_TRUE(tokens.at(1).empty());
	}
	
	{
		std::string s = ",b";
		auto tokens = split(s,',','(',')');
		ASSERT_TRUE(tokens.size()==2u);
		ASSERT_TRUE(tokens.at(0).empty());
		ASSERT_TRUE(tokens.at(1) == "b");
	}
	
	{
		std::string s = "a,b";
		auto tokens = split(s,',','(',')');
		ASSERT_TRUE(tokens.size()==2u);
		ASSERT_TRUE(tokens.at(0) == "a");
		ASSERT_TRUE(tokens.at(1) == "b");
	}
	
	{
		std::string s = "(a,b)";
		auto tokens = split(s,',','(',')');
		ASSERT_TRUE(tokens.size()==1u);
		ASSERT_TRUE(tokens.at(0) == "(a,b)");		
	}
	
	{
		std::string s = "(a,b),";
		auto tokens = split(s,',','(',')');
		ASSERT_TRUE(tokens.size()==2u);
		ASSERT_TRUE(tokens.at(0) == "(a,b)");		
		ASSERT_TRUE(tokens.at(1).empty());
	}
	
	{
		std::string s = ",(a,b)";
		auto tokens = split(s,',','(',')');
		ASSERT_TRUE(tokens.size()==2u);
		ASSERT_TRUE(tokens.at(0).empty());		
		ASSERT_TRUE(tokens.at(1) == "(a,b)");
	}
	
	{
		std::string s = "(a,b),(c,d)";
		auto tokens = split(s,',','(',')');
		ASSERT_TRUE(tokens.size()==2u);
		ASSERT_TRUE(tokens.at(0) == "(a,b)");		
		ASSERT_TRUE(tokens.at(1) == "(c,d)");
	}
	
	{
		std::string s = "()";
		auto tokens = split(s,',','(',')');
		ASSERT_TRUE(tokens.size()==1u);
		ASSERT_TRUE(tokens.at(0)=="()");				
	}
	
	{
		std::string s = ",(),";
		auto tokens = split(s,',','(',')');
		ASSERT_TRUE(tokens.size()==3u);
		ASSERT_TRUE(tokens.at(0).empty());				
		ASSERT_TRUE(tokens.at(1)=="()");				
		ASSERT_TRUE(tokens.at(2).empty());				
	}
	
	{ // not recursiv
		std::string s = "(a,(b,c),d)";
		auto tokens = split(s,',','(',')');
		ASSERT_TRUE(tokens.size()==2u);
		ASSERT_TRUE(tokens.at(0) == "(a,(b,c)");		
		ASSERT_TRUE(tokens.at(1) == "d)");
	}
}



TEST(quote,test)
{
	{
		auto s = quote("");
		ASSERT_TRUE(s=="\"\"");
	}
	{
		auto s = quote("a");
		ASSERT_TRUE(s=="\"a\"");
	}
}



TEST(differences,test)
{
	{
		std::array<double,1> values {1};
		std::array<double,0> expected {};
		auto result = differences(values);
		ASSERT_TRUE(result==expected);
	}

	{
		std::array<double,2> values {1,1};
		std::array<double,1> expected {0};
		auto result = differences(values);
		ASSERT_TRUE(result==expected);
	}

	{
		std::array<double,4> values {1,2,4,5};
		std::array<double,3> expected {1,2,1};
		auto result = differences(values);
		ASSERT_TRUE(result==expected);
	}
}


TEST(starts_with,test)
{
	{
		ASSERT_TRUE(starts_with("",""));
		ASSERT_FALSE(starts_with("a",""));
		ASSERT_FALSE(starts_with("","a"));
		ASSERT_TRUE(starts_with("a","a"));
		ASSERT_FALSE(starts_with("b","a"));
		ASSERT_TRUE(starts_with("ab","a"));
		ASSERT_FALSE(starts_with("ba","a"));
		ASSERT_FALSE(starts_with("a","ab"));
		ASSERT_FALSE(starts_with("a","ba"));
	}
}


TEST(ends_with,test)
{
	{
		ASSERT_TRUE(ends_with("",""));
		ASSERT_FALSE(ends_with("","a"));
		ASSERT_FALSE(ends_with("a",""));
		ASSERT_TRUE(ends_with("a","a"));
		ASSERT_FALSE(ends_with("a","b"));
		ASSERT_FALSE(ends_with("a","ba"));
		ASSERT_FALSE(ends_with("a","ab"));
		ASSERT_TRUE(ends_with("ba","a"));
		ASSERT_FALSE(ends_with("ab","a"));
	}
}


TEST(trim,test)
{
	{
		ASSERT_TRUE(trim("") == "");
		ASSERT_TRUE(trim("\t") == "");
		ASSERT_TRUE(trim("\n") == "\n");
		ASSERT_TRUE(trim("a") == "a");
		ASSERT_TRUE(trim("a ") == "a");
		ASSERT_TRUE(trim("a  ") == "a");
		ASSERT_TRUE(trim(" a") == "a");
		ASSERT_TRUE(trim("  a") == "a");
		ASSERT_TRUE(trim(" a a ") == "a a");
	}
}






TEST(parse_double,test)
{
	try {
		parse_double("");
		FAIL();
	}
	catch(std::invalid_argument&) {
	}

	try {
		parse_double(" ");
		FAIL();
	}
	catch(std::invalid_argument&) {
	}

	try {
		parse_double("-");
		FAIL();
	}
	catch(std::invalid_argument&) {
	}

	try {
		parse_double("+");
		FAIL();
	}
	catch(std::invalid_argument&) {
	}

	try {
		parse_double(".");
		FAIL();
	}
	catch(std::invalid_argument&) {
	}

	try {
		parse_double(",");
		FAIL();
	}
	catch(std::invalid_argument&) {
	}

	try {
		parse_double("0,1");
		FAIL();
	}
	catch(std::invalid_argument&) {
	}

	try {
		parse_double("0.1 ");
		FAIL();
	}
	catch(std::invalid_argument&) {
	}




	{
		ASSERT_TRUE(parse_double("1")==1);
		ASSERT_TRUE(parse_double("0.1")==0.1);
		ASSERT_TRUE(parse_double(" 0.1")==0.1);
		ASSERT_TRUE(parse_double("-0.1")==-0.1);
		ASSERT_TRUE(parse_double("1.")==1);
		ASSERT_TRUE(parse_double(".1")==0.1);
	}
}



TEST(max_abs_value,test)
{
	{
		std::array<double,0> a;
		auto m = max_abs_value(a);
		ASSERT_TRUE(m==0);
	}

	{
		std::array<double,1> a {-1};
		auto m = max_abs_value(a);
		ASSERT_TRUE(m==1);
	}


	{
		std::array<double,5> a {-1,2,3,-2,0};
		auto m = max_abs_value(a);
		ASSERT_TRUE(m==3);
	}

	{
		std::array<double,5> a {-1,2,-3,-2,0};
		auto m = max_abs_value(a);
		ASSERT_TRUE(m==3);
	}
}




TEST(add_array,test)
{
	{
		std::array<double,3> a {};
		std::array<double,3> b {1,2,3};
		add_array(a,b);
		ASSERT_TRUE(a==b);
	}


	{
		std::array<double,3> a {3,2,1};
		std::array<double,3> b {-1,-2,-3};
		std::array<double,3> expected {2,0,-2};
		add_array(a,b);
		ASSERT_TRUE(a==expected);
	}
}




TEST(sub_array,test)
{
	{
		std::array<double,3> a {};
		std::array<double,3> b {1,2,3};
		std::array<double,3> expected {-1,-2,-3};
		sub_array(a,b);
		ASSERT_TRUE(a==expected);
	}


	{
		std::array<double,3> a {3,2,1};
		std::array<double,3> b {-1,-2,-3};
		std::array<double,3> expected {4,4,4};
		sub_array(a,b);
		ASSERT_TRUE(a==expected);
	}
}




TEST(add_sqr_array,test)
{
	{
		std::array<double,3> a {};
		std::array<double,3> b {1,2,3};
		std::array<double,3> expected {1,4,9};
		add_sqr_array(a,b);
		ASSERT_TRUE(a==expected);
	}


	{
		std::array<double,3> a {3,2,1};
		std::array<double,3> b {-1,-2,-3};
		std::array<double,3> expected {3+1,2+4,1+9};
		add_sqr_array(a,b);
		ASSERT_TRUE(a==expected);
	}
}



TEST(for_each_of,test)
{
	{
		std::array<double,3> a {1.1, 2.2, 3.3};
		std::vector<int> v;

		int n=0;
		double sum=0;
		for_each_of(a,v,[&n,&sum](double ai,int vi){
			sum += ai+vi;
			n++;
		});
		ASSERT_TRUE(n==0 && sum==0);
	}

	{
		std::array<double,3> a {1.1, 2.2, 3.3};
		std::vector<int> v {1,2,3,4};

		int n=0;
		double sum=0;
		for_each_of(a,v,[&n,&sum](double ai,int vi){
			sum += ai+vi;
			n++;
		});
		ASSERT_TRUE(n==3);
		ASSERT_APPROX(sum,12.6,0.001);
	}
}



TEST(time_of_day,test)
{
	{
		ASSERT_TRUE(time_of_day(7,19,6.9)==TimeOfDay::Night);
		ASSERT_TRUE(time_of_day(7,19,7)==TimeOfDay::Sunrise);
		ASSERT_TRUE(time_of_day(7,19,7+1)==TimeOfDay::Day);
		ASSERT_TRUE(time_of_day(7,19,19)==TimeOfDay::Night);
		ASSERT_TRUE(time_of_day(7,19,19-1)==TimeOfDay::Sunset);
		ASSERT_TRUE(time_of_day(7,19,19-1.1)==TimeOfDay::Day);
	}
}



TEST(create_random_seed,test)
{
	{
		auto s1 = create_random_seed();
		auto s2 = create_random_seed();
		ASSERT_TRUE(s1!=s2);
	}
}



TEST(julian_date,exceptions)
{
	{
		Calendar d;
		d.year=2019;
		d.month=0;
		d.day=1;
		try {
			julian_date(d);
			FAIL();
		}
		catch(std::invalid_argument& ) {
		}
	}

	{
		Calendar d;
		d.year=2019;
		d.month=13;
		d.day=1;
		try {
			julian_date(d);
			FAIL();
		}
		catch(std::invalid_argument& ) {
		}
	}

	{
		Calendar d;
		d.year=2019;
		d.month=1;
		d.day=0;
		try {
			julian_date(d);
			FAIL();
		}
		catch(std::invalid_argument& ) {
		}
	}
}



TEST(julian_date,test)
{
	{
		Calendar d {1970,1,1};
		auto jd = julian_date(d);
		ASSERT_APPROX(jd, 2440587.5,0.001);
	}

	{
		Calendar d {-4712,1,1 + (12./24)};
		auto jd = julian_date(d,false);
		ASSERT_APPROX(jd, 0,0.001);
	}

	{
		Calendar d {-4713,11,24.5};
		auto jd = julian_date(d);
		ASSERT_APPROX(jd, 0,0.001);
	}
}



TEST(julian_date,wikipedia)
{
	{
		// https://de.wikipedia.org/wiki/Julianisches_Datum
		Calendar d {-668,5,27 + (1 + 59/60.0)/24.0};
		auto jd = julian_date(d,false);
		
		ASSERT_APPROX(jd,1477217.583,0.01); 
	}

	{
		// https://de.wikipedia.org/wiki/Julianisches_Datum
		Calendar d {1,1,1 + (0)};
		auto jd = julian_date(d,false);

		ASSERT_APPROX(jd,1721423.500, 0.01); 
	}



	{
		// https://en.wikipedia.org/wiki/julian_date#cite_ref-8
		// US Naval Observatory 2005:

		Calendar d {2013,1,1 + (0 + 30/60.0)/24.0};
		auto jd = julian_date(d);

		ASSERT_APPROX(jd,2456293.520833, 0.00001);
	}
}



TEST(parse_ulong,test)
{
	{
		try {
			parse_ulong("");
			FAIL();
		}
		catch(std::invalid_argument&) {			
		}
	}
	
	{
		try {
			parse_ulong("a");
			FAIL();
		}
		catch(std::invalid_argument&) {			
		}
	}
	
	{
		try {
			auto ul = parse_ulong("-0"); // will be parsed			
			//FAIL();
			ASSERT_TRUE(ul==0);		
		}
		catch(std::invalid_argument&) {			
		}
	}
	
	{
		try {
			parse_ulong("-1");
			FAIL();
		}
		catch(std::invalid_argument&) {			
		}
	}
	
	{		
		auto ul = parse_ulong("0");
		ASSERT_TRUE(ul==0);		
	}
	
	{		
		auto ul = parse_ulong(" 0");
		ASSERT_TRUE(ul==0);		
	}	
	
	{		
		auto ul = parse_ulong("13");
		ASSERT_TRUE(ul==13);		
	}
	
	{		
		auto ul = parse_ulong("4100000000");
		ASSERT_TRUE(ul==4100000000u);		
	}
	
	
}


int main(int ,char** )
{
	RUN_ALL_TESTS();

	return EXIT_SUCCESS;
}
