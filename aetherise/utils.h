/**
 * \file
 *
 * \~german
 * @brief Allgemeines und Nützliches
 *
 * \~english
 * @brief Utilities and common stuff
 *
 */

#ifndef UTILS_H
#define UTILS_H

#include "stdx.h"
#include "mathematics.h"

#include <string>
#include <vector>
#include <memory>
#include <cmath>
#include <algorithm>
#include <ostream>
#include <sstream>
#include <unordered_set>



#define CSV_DELIMITER ';'
#define DATE_TIME_ISO "%Y-%m-%d %H:%M:%S"

#define SSSTR(p) (static_cast<std::ostringstream&>(std::ostringstream() << p)).str()
#define WSSTR(p) (static_cast<std::wostringstream&>(std::wostringstream() << p)).str()




namespace aether {



constexpr const char* yesno(bool b)
{
	return b ? "yes" : "no";	
}



extern bool locale_german;

struct _StoreCursor{}; // helper
struct _RestoreCursor{}; // helper

std::ostream& operator << (std::ostream& os,const _StoreCursor&);
std::ostream& operator << (std::ostream& os,const _RestoreCursor&);

const _StoreCursor store_cursor;
const _RestoreCursor restore_cursor;



template<typename T>
std::ostream& operator <<(std::ostream& os,const optional<T>& o)
{
	if (o.has_value())
		os << *o;
	else
		os << ""; // needed for std::setw() to work correctly
	return os;
}



/**
 * \~german
 * Versionsnummer
 *
 * \~english
 * Version number
 */
struct Version
{
	short major,minor,patch;
};

std::ostream& operator <<(std::ostream& os,const Version& version);



/**
 * \~german
 * Fehler bei der Textverarbeitung
 *
 * \~english
 * Error while processing text
 */
struct ParseException
{
	std::string message;
	int pos;
};



struct ExitException
{
	std::string message;
};



/**
 * \~german
 * Parameter der Theorie
 *
 * Der Vektor (v,a,d) gibt die Geschwindigkeit und Richtung im Äther an.
 *
 * \~english
 * Parameters of the theory
 *
 * The vector (v,a,d) specifies the speed and direction of movement in the aether.
 *
 */
struct TheoryParameters
{
	double v; ///< velocity (m/s)
	double a; ///< right ascension (rad)
	double d; ///< declination (rad)
};



// Direction of movement according to the CMB dipole measured by WMAP.
// (galactic coordinates: l=263.99 b=48.26)
const double CMB_DIPOLE_V = 369000; // (m/s) WMAP velocity value of CMB dipole
const double CMB_DIPOLE_RA = h_to_rad(11.195); // 11h 11m
const double CMB_DIPOLE_DE = rad(-6.93); // -6° 55'

const TheoryParameters CMB_dipole {CMB_DIPOLE_V, CMB_DIPOLE_RA, CMB_DIPOLE_DE};

const double MtWilson_Latitude = rad(34.225);
const double MtWilson_Longitude = rad(-118.057);
const double MtWilson_Altitude = 1700; // (m)

const double Millers_Interferometer_Wave_Length = 570e-9; // (m) acetylene lamp
const double Millers_Interferometer_Arm_Length = 32.03; // (m) 


struct IntegerInterval
{
	int from;
	int to;
};


/**
 * \~german
 * Mögliche Arten der Datenverarbeitung
 *
 * \~english
 * Possible kinds of data processing
 *
 */
enum class Action
{
	Filename,Header,Raw,RawReduced,Reduce,Test,Spectrum,RawSpectrum
};



/**
 * \~german
 * Einstellungen zur Datenverarbeitung
 *
 * \~english
 * Options for data processing
 */
struct Options
{
	enum class ReductionMethod {
		Miller,DFT,Roberts2006
	};

	enum class OutputFormat {
		CSV,Gnuplot
	};

	enum class Theory {
		Classic,Aether,Relativity
	};

	enum class AggregationMethod {
		List,Mean,Sidereal,Diff,Params,ModelChi,Fit,Test,Signals
	};

	enum class Minimizer {
		Grad, Minuit2
	};

	Theory theory = Theory::Aether;
	AggregationMethod aggregation_method = AggregationMethod::Mean;
	ReductionMethod reduction_method = ReductionMethod::Miller;
	Minimizer minimizer = Minimizer::Grad;
	OutputFormat output_format = OutputFormat::Gnuplot;
	std::string data_filename;
	std::vector<double> data1;
	std::vector<double> data2;
	std::vector<double> data3;
	std::vector<double> data4;
	double delta_chi_squared = 1.0;	
	double chi_squared_scale = 1.0;
	double signals_dTD = 0.1;
	double signals_ddT = 0.1;
	double signals_dt = 1.0;
	optional<TheoryParameters> start_params;
	TheoryParameters theory_params {CMB_dipole};
	std::unordered_set<int> disabled_signals;
	double index_of_refraction = NAN; // not set
	double latitude = MtWilson_Latitude; 
	double longitude = MtWilson_Longitude; 
	double altitude = MtWilson_Altitude;
	double arm_length = Millers_Interferometer_Arm_Length;	
	unsigned int sim_seed = 0;
	bool invert_data = false;
	bool invert_theory = false;
	bool invert_model = false;
	bool invert_temp = false;
	bool add = false;
	bool subtract = false;
	bool validate = false;
	bool subtract_data = false;
	bool subtract_model = false;	
	bool single = false;	
	bool stats = false;
	bool model = false;
	bool enable_desk = true;
	bool enable_temp = true;
	bool enable_earth = true;
	bool output_data = true;
	bool output_theory = true;
	bool epoch_params = false;
	bool fit_amplitude = false;	
	bool fit_sine = false;	
	bool loocv = false;	
	bool ignore_reverse = false;
	bool ignore_cancel = false;
	bool ignore_reverse_sheet = false;
	bool ignore_invert = false;
	bool ignore_bad = false;
	bool ignore_reverse_disabled = false;
	bool ignore_cancel_disabled = false;
	bool ignore_sw = false;
	bool ignore_inverted_theory = false;
	bool contour = false;
	bool day_and_night = false;
	bool low_sun = false;
	bool simulation = false;
	bool sim_simple = false;
	bool sim_sys = false;
	bool residuals = false;
	bool subtract_first_half_turn = true;	
};




/**
 * \~german Tageszeit
 *
 * Alles was nicht Nacht ist, ist Tag.
 *
 * \~english
 * Everything that is not night, is day.
 */
enum class TimeOfDay
{
	Night,Sunrise,Day,Sunset
};



/**
 * \~german
 * Tageszeit aus dem Zeitpunkt von Sonnenaufgang und Sonnenuntergang ermitteln.
 * Sonnenaufgang und Sonnenuntergang dauern 1 Stunde.
 * Alle Parameter erwarten Argumente im Intervall [0, 24).
 *
 * \~english
 * Evalute time of day for given times of sunrise and sunset.
 * Sunrise and sunset last 1 hour.
 * All parameters expect arguments inside interval [0, 24).
 *
 * \~
 * @param sunrise
 * @param sunset
 * @param time
 * @return
 */
TimeOfDay time_of_day(double sunrise,double sunset,double time);



struct Calendar
{
	int year;
	short month;
	double day;
};




template<class C1,class C2,typename F>
void for_each_of(C1&& c1, C2&& c2,F f)
{
	auto i1 = std::begin(std::forward<C1>(c1));
	auto i2 = std::begin(std::forward<C2>(c2));
	for (;i1!=std::end(std::forward<C1>(c1)) && i2!=std::end(std::forward<C2>(c2)); ++i1,++i2) {
		f(*i1,*i2);
	}
}



/**
 * \~german
 * Den größten absoluten Wert finden.
 *
 * \~english
 * Find the greatest absolute value.
 *
 * \~
 * @return
 */
template<typename T,size_t N>
double max_abs_value(const std::array<T,N>& values)
{
	T max_value {0};
	for (const T& value : values) {
		max_value = std::max(max_value,std::abs(value));
	}
	return max_value;
}



/**
 * \~german Differenzen
 *
 * Die Differenzen benachbarter Elemente.
 *
 * \~english
 * The differences of adjacent elements.
 *
 * \~
 * @param values
 * @return
 */
template<typename T,size_t N>
std::array<T,N-1> differences(const std::array<T,N>& values)
{
	std::array<T,N-1> ds;
	auto it = values.begin();
	for (T& d : ds) {
		auto next = it+1;
		d = (*next - *it);
		it = next;
	}
	return ds;
}




/**
 * a <- a+b
 */
template<typename T,size_t N>
void add_array(std::array<T,N>& a,const std::array<T,N>& b)
{
	std::transform(a.begin(),a.end(),b.begin(),a.begin(),std::plus<T>());
}

/**
 * a <- a-b
 */
template<typename T,size_t N>
void sub_array(std::array<T,N>& a,const std::array<T,N>& b)
{
	std::transform(a.begin(),a.end(),b.begin(),a.begin(),std::minus<T>());
}

/**
 * a <- a+b²
 */
template<typename T,size_t N>
void add_sqr_array(std::array<T,N>& a,const std::array<T,N>& b)
{
	std::transform(a.begin(),a.end(),b.begin(),a.begin(),[](T ai,T bi){
		return ai + bi*bi;
	});
}



/**
 * \~german
 * Eine Sequenz in Teile zerlegen.
 * Funktor f wird mindestens einmal aufgerufen.
 *
 * \~english
 * Split a sequence.
 * Functor f is called at least once.
 */
template<typename Iter,typename T,typename Functor>
void split(Iter begin,Iter end,T delimiter,Functor f)
{
	using I = typename std::iterator_traits<Iter>::value_type;
	static_assert(std::is_same<I,T>::value,"must be same types");
	
	Iter it = begin;
	while (it != end) {
		if (*it == delimiter) {
			f(begin,it);
			begin = ++it;
		}
		else
			++it;
	}

	f(begin,it);
}


/**
 * \~german
 * Eine Sequenz in Teile zerlegen.
 * Trenner zwischen den gegebenen Klammerwerten haben keine Wirkung.
 * Funktor f wird mindestens einmal aufgerufen.
 *
 * \~english
 * Split a sequence.
 * Delimiters between the given bracket values have no effect.
 * Functor f is called at least once.
 */
template<typename Iter,typename T,typename Functor>
void split(Iter begin,Iter end,T delimiter,T left_bracket,T right_bracket,Functor f)
{
	using I = typename std::iterator_traits<Iter>::value_type;
	static_assert(std::is_same<I,T>::value,"must be same types");
	
	bool in_text = false;
	Iter it = begin;
	while (it != end) {
		if (!in_text && *it == delimiter) {
			f(begin,it);
			begin = ++it;
		}
		else if (!in_text && *it == left_bracket) {
			in_text = !in_text;
			++it;
		}
		else if (in_text && *it == right_bracket) {
			in_text = !in_text;
			++it;
		}
		else
			++it;
	}

	f(begin,it);
}



/**
 * \~german
 * Zeichenkette in Teile zerlegen.
 * Trennzeichen zwischen doppelten Anführungsstrichen ("") haben keine Wirkung.
 *
 * \~english
 * Split char squence into parts.
 * Delimiters between double quotes ("") have no effect.
 *
 * \~
 * @param s
 * @param delimiter
 * @return tokens
 */
std::vector<std::string> split(const std::string& s,char delimiter);



/**
 * \~german
 * Eine Sequenz in Teile zerlegen.
 * Trennzeichen zwischen den gegebenen Klammerzeichen haben keine Wirkung.
 * 
 * \~english
 * Split a sequence.
 * Delimiters between the given bracket chars have no effect.
 * 
 * \~
 * @param s
 * @param delimiter
 * @param left_bracket
 * @param right_bracket
 * @return tokens
 */
std::vector<std::string>
split(const std::string& s,char delimiter,char left_bracket,char right_bracket);



/**
 * \~german
 * Zeichen vom Anfang und Ende der gegebenen Zeichenkette entfernen.
 *
 * \~english
 * Remove chars from beginning and end of a string.
 *
 * \~
 * @param s
 * @param chars
 * @return
 */
template<typename T>
std::basic_string<T>
trim(const std::basic_string<T>& s,const std::basic_string<T>& chars)
{
	auto p0 = s.find_first_not_of(chars);
	if (p0 == std::basic_string<T>::npos)
		return {};
	auto p1 = s.find_last_not_of(chars);
	return s.substr(p0,p1-p0+1);
}



/**
 * \~german
 * Leerzeichen und Tabulator vom Anfang und Ende entfernen.
 *
 * \~english
 * Remove spaces and tabs from beginning and end.
 *
 * \~
 * @param s
 * @return
 */
std::string trim(const std::string& s);



/**
 * \~german
 * Prüfen ob ein String mit der gegebenen Zeichenkette anfängt.
 *
 * \~english
 * Check if a string starts with the given char sequence.
 *
 * \~
 * @param str
 * @param seq
 * @return
 */
bool starts_with(const std::string& str,const std::string& seq);



/**
 * \~german
 * Prüfen ob ein String mit der gegebenen Zeichenkette endet.
 *
 * \~english
 * Check if a string ends with the given char sequence.
 *
 * \~
 * @param str
 * @param seq
 * @return
 */
bool ends_with(const std::string& str,const std::string& seq);






template<typename Iter,typename S,typename Functor>
void output_separated(std::ostream& os,Iter it,Iter end,S separator,Functor get)
{
	if (it == end)
		return;
	os << get(*it);
	for(++it; it!=end; ++it) {
		os << separator << get(*it);
	}
}

template<typename Iter,typename S>
void output_separated(std::ostream& os,Iter it,Iter end,S separator)
{
	using T = typename std::iterator_traits<Iter>::value_type;

	output_separated(os,it,end,separator,[](T& v){
		return v;
	});
}


template<class C,typename S,typename Functor>
void output_separated(std::ostream& os,const C& c,S separator,Functor get)
{	
	output_separated(os,std::begin(c),std::end(c),separator,get);
}


template<class C,typename S>
void output_separated(std::ostream& os,const C& c,S separator)
{
	using T = decltype(*std::begin(c));

	output_separated(os,c,separator,[](T& v){
		return v;
	});
}




template<typename S,class Arg>
void __output_args_separated(std::ostream& os,S ,Arg&& arg)
{
	os << std::forward<Arg>(arg);
}

template<typename S,class Arg,class... Args>
void __output_args_separated(std::ostream& os,S separator,Arg&& arg,Args&&... args)
{
	os << std::forward<Arg>(arg) << separator;
	__output_args_separated(os,separator,std::forward<Args>(args)...);
}

/**
 * @param args
 */
template<typename S,class... Args>
std::ostream& output_args_separated(std::ostream& os,S separator,Args&&... args)
{
	__output_args_separated(os,separator,std::forward<Args>(args)...);
	return os;
}



/**
 * \~german
 * Verschiedene gnuplot "data sets" in einer Datei trennen.
 * Wird z.B. für die 3D Ansicht benötigt, damit keine Verbindunglinien zwischen
 * einzelnen Kurven gezeichnet werden.
 *
 * \~english
 * Separate different data sets in a gnuplot data file.
 * Used for the 3D view, to not connect the separate curves.
 *
 * \~
 * @param os
 */
inline void write_gnuplot_data_set_separator(std::ostream& os)
{
	os << "\n\n";
}




/**
 * \~german
 * Julianisches Datum (JD) aus gegebenen Kalenderdatum berechnen.
 *
 * \~english
 * Julian date (JD) calculated from given calendar.
 *
 * \~
 * @param cal
 * @param gregorian Gregorian calendar? Set false to use Julian calendar.
 * @return JD
 */
double julian_date(Calendar cal,bool gregorian=true);






/**
 * \~german
 * Zeichenkette nach Gleitkommatyp double wandeln.
 * Es wird die englische Schreibweise erwartet, also Punkt statt Komma.
 *
 * \~english
 * Converts a string, containing a decimal, into a double.
 *
 * \~
 * @param str
 * @return
 * @throws invalid_argument
 */
double parse_double(const std::string& str);



/**
 * \~german
 * Zeichenkette nach Ganzzahltyp int wandeln.
 *
 * \~english
 * Converts a string, containing an integer, into an int.
 *
 * \~
 * @param str
 * @return
 * @throws invalid_argument
 */
int parse_int(const std::string& str);



/**
 * \~german
 * Zeichenkette nach Ganzzahltyp unsigned long wandeln.
 * 
 * \~english
 * Converts a string, containing an integer, into an unsigned long.
 *
 * \~
 * @param str
 * @return 
 */
unsigned long parse_ulong(const std::string& str);



/**
 * \~german
 * Zitat/zitieren
 *
 * Gegebene Zeichenkette in doppelte Anführungstriche setzten.
 *
 * \~english
 * Double quote the given string.
 *
 * \~
 * @param s
 * @return
 */
std::string quote(const std::string& s);



/**
 * \~german
 * Zufälligen Startwert für Zufallszahlengeneratoren erzeugen.
 * Auf std::random_device kann man sich nicht verlassen!
 *
 * \~english
 * Create a random seed for random number engines.
 * std::random_device is not reliable!
 *
 * \~
 * @return
 */
unsigned int create_random_seed();

}//aether

extern const aether::Version version;


#endif // UTILS_H
