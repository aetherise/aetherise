/**
 * \file
 *
 * \~german
 * @brief Einfaches Testgerüst
 *
 * \~english
 * @brief Simple test framework
 *
 */


#ifndef TEST_H
#define TEST_H


#include <cstdlib>
#include <iostream>
#include <vector>
#include <functional>

namespace test {

//using _Test = std::function<void()>;
struct _Test {
	const char* name1;
	const char* name2;
	std::function<void()> f;

	void operator ()(void) const {
		std::cout << "TEST: " << name1 << ": " << name2 << "\n";
		f();
	}

};

extern std::vector<_Test> _tests;

class _add_test {
public:
	_add_test(const _Test& t) {
		_tests.push_back(t);
	}
};

bool approximate(double a, double b,double f = 0.000001);

}//test



#define TEST(name1,name2) \
	void name1##__##name2();\
	static test::_add_test _test_##name1##__##name2({#name1,#name2,name1##__##name2});\
	void name1##__##name2()

#define RUN_ALL_TESTS(null) for (const test::_Test& test : test::_tests) test();



#define ASSERT_TRUE(condition)\
	if (!(condition)) {\
		std::cerr << "Test " << __func__ << "() failed in " << __FILE__ << ":" << __LINE__ << "\n";\
		std::cerr << "Assertion '" << #condition << "' failed\n";\
		std::exit(EXIT_FAILURE);\
	}

#define ASSERT_FALSE(condition)\
	ASSERT_TRUE(!(condition));


#define ASSERT_APPROX(val1,val2,precision)\
	if (!test::approximate(val1,val2,precision)) {\
		std::cerr << "Test " << __func__ << "() failed in " << __FILE__ << ":" << __LINE__ << "\n";\
		std::cerr << val1 << " is not approximate " << val2 << " using precision " << precision << "\n";\
		std::exit(EXIT_FAILURE);\
	}


#define FAIL(msg)\
	std::cerr << "Test " << __func__ << "() failed in " << __FILE__ << ":" << __LINE__ << "\n";\
	std::cerr << "" msg << "\n";\
	std::exit(EXIT_FAILURE);


#endif // TEST_H
