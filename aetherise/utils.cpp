
#include "utils.h"

#ifdef AETHER_WINDOWS
#include <windows.h>
#endif

#include <stdexcept>
#include <cmath>
#include <iomanip>
#include <ctime>
#include <ostream>
#include <sstream>
#include <iostream>
#include <random>


namespace aether {

bool locale_german {false};

std::ostream& operator <<(std::ostream& os,const Version& version)
{
	std::ostringstream oss;
	oss << version.major << '.' << version.minor << '.' << version.patch;
	return os << oss.str();
}




#ifdef AETHER_WINDOWS


const int StoreRestoreCursor_id = std::ios_base::xalloc();

std::ostream& operator << (std::ostream& os,const _StoreCursor&) {
	auto H = &os==&std::cerr ? STD_ERROR_HANDLE : STD_OUTPUT_HANDLE;
	auto h = ::GetStdHandle(H);
	if (h!=INVALID_HANDLE_VALUE) {
		CONSOLE_SCREEN_BUFFER_INFO csbi;
		if(::GetConsoleScreenBufferInfo(h,&csbi)) {
			static_assert(sizeof(COORD::X)<=2,"16 bit SHORT expected");
			static_assert(sizeof(COORD::Y)<=2,"16 bit SHORT expected");
			static_assert(sizeof(int)>=4,"32 bit int needed");
			os.iword(StoreRestoreCursor_id) = (csbi.dwCursorPosition.X << 16) | csbi.dwCursorPosition.Y;
		}
	}

	return os;
}

std::ostream& operator << (std::ostream& os,const _RestoreCursor&) {
	auto H = &os==&std::cerr ? STD_ERROR_HANDLE : STD_OUTPUT_HANDLE;
	auto h = ::GetStdHandle(H);
	if (h!=INVALID_HANDLE_VALUE) {
		int iw = os.iword(StoreRestoreCursor_id);
		COORD coord;
		coord.X = iw >> 16;
		coord.Y = iw & 65535;
		::SetConsoleCursorPosition(h,coord);
	}

	return os;
}


#else


std::ostream& operator << (std::ostream& os,const _StoreCursor&) {
	return os << "\x1b[s"; // ANSI escape sequence
}

std::ostream& operator << (std::ostream& os,const _RestoreCursor&) {
	return os << "\x1b[u"; // ANSI escape sequence
}

#endif




TimeOfDay time_of_day(double sunrise,double sunset,double time)
{
	if (time < sunrise)
		return TimeOfDay::Night;
	if (time >= sunset)
		return TimeOfDay::Night;

	if (time < sunrise + 1)
		return TimeOfDay::Sunrise;
	if (time >= sunset - 1)
		return TimeOfDay::Sunset;

	return TimeOfDay::Day;
}




std::vector<std::string> split(const std::string& s,char delimiter)
{
	using CIter = std::string::const_iterator;
	std::vector<std::string> tokens;

	split(s.begin(),s.end(),delimiter,'"','"',[&tokens](CIter from,CIter to) {
		tokens.emplace_back(from,to);
	});

	return tokens;
}


std::vector<std::string>
split(const std::string& s,char delimiter,char left_quote,char right_quote)
{
	using CIter = std::string::const_iterator;
	std::vector<std::string> tokens;

	split(s.begin(),s.end(),delimiter,left_quote,right_quote,[&tokens](CIter from,CIter to) {
		tokens.emplace_back(from,to);
	});

	return tokens;
}

std::string trim(const std::string& s)
{
	return trim(s,std::string(" \t"));
}




bool starts_with(const std::string& str,const std::string& chars)
{
	if (!str.empty() && chars.empty())
		return false;
	return std::equal(chars.begin(),chars.end(),str.begin());
}



bool ends_with(const std::string& str,const std::string& chars)
{
	if (!str.empty() && chars.empty())
		return false;
	return std::equal(chars.rbegin(),chars.rend(),str.rbegin());
}



double julian_date(Calendar cal,bool gregorian)
{
	if (cal.month<1 || cal.month>12)
		throw std::invalid_argument("calendar month");
	if (cal.day<1)
		throw std::invalid_argument("calendar day");

	// Peter Baum, Date Algorithms
	// https://www.researchgate.net/publication/316558298_Date_Algorithms

	if (cal.month < 3) {
		cal.month += 12;
		cal.year -= 1;
	}

	double j = 365.0*cal.year + std::floor(cal.year/4.0);

	double t;
	if (gregorian) {
		j += std::floor(cal.year/400.0) - std::floor(cal.year/100.0);
		t = 1721118.5;
	}
	else
		t = 1721116.5;

	int m = (979*cal.month-2918)/32; // integer division
	return cal.day + m + j + t;
}




double parse_double(const std::string& str)
{
	try {
		size_t n = 0;
		double d = std::stod(str,&n);
		if (n != str.size())
			throw std::invalid_argument("not a decimal");
		return d;
	}
	catch(std::exception& ) {
		throw std::invalid_argument("not a decimal");
	}
}



int parse_int(const std::string& str)
{
	try {
		size_t n = 0;
		int i = std::stoi(str,&n);
		if (n != str.size())
			throw std::invalid_argument("not an integer");
		return i;
	}
	catch(std::exception& ) {
		throw std::invalid_argument("not an integer");
	}
}




unsigned long parse_ulong(const std::string& str)
{
	try {
		size_t n = 0;
		unsigned long ul = std::stoul(str,&n);
		if (n != str.size())
			throw std::invalid_argument("not an positive integer");
		if (str.find("-")!=std::string::npos) // stoul parses negative values
			throw std::invalid_argument("not an positive integer");
		return ul;
	}
	catch(std::exception& ) {
		throw std::invalid_argument("not an positive integer");
	}
}



std::string quote(const std::string& s)
{
	std::string q;
	q.reserve(2+s.size());
	q.append("\"").append(s).append("\"");
	return q;
}



unsigned int create_random_seed()
{	
	static unsigned random_counter; // not initialized!
	
	// random_device is not reliable!
	// It can behave constant or deterministic on some compilers and systems!
	std::random_device rdev;
	auto t = std::time(nullptr);	
	std::seed_seq seq {
				unsigned(rdev()),unsigned(t),unsigned(std::clock()),unsigned(random_counter++)				
	};
	std::array<unsigned,1> numbers;
	seq.generate(numbers.begin(),numbers.end());
	return numbers[0];
}

}//aether
