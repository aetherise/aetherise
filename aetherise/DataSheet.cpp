
#include "DataSheet.h"
#include "utils.h"
#include "mathematics.h"
#include "astro.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <limits>
#include <algorithm>
#include <cmath>

#define HEADER_LINES 4


#define USUAL_NUMBER_OF_TURNS 20
#define MAX_DURATION 45
#define MIN_DURATION 10
#define MAX_TEMP_DIFF 2

namespace aether {



bool operator < (const DataSheet::Date& a, const DataSheet::Date& b)
{
	if (a.year == b.year) {
		if (a.month == b.month) {
			return a.day < b.day;
		}
		return a.month < b.month;
	}

	return a.year < b.year;
}



bool operator == (const DataSheet::Date& a, const DataSheet::Date& b)
{
	return a.year==b.year && a.month==b.month && a.day==b.day;
}



bool operator < (const DataSheet::Time& a,const DataSheet::Time& b)
{
	if (a.hour == b.hour)
		return a.minute < b.minute;
	return a.hour < b.hour;
}


bool operator == (const DataSheet::Time& a,const DataSheet::Time& b)
{
	return a.hour==b.hour && a.minute==b.minute;
}


bool operator <(const DataSheet& a, const DataSheet& b)
{
	if (a.date == b.date)
		return a.no < b.no;
	return a.date < b.date;
}



std::ostream& operator<<(std::ostream& os,const DataSheet::Date& date)
{
	std::ostringstream ss; // setw should work
	ss << date.year << "-";
	if(date.month < 10)
		ss << "0";
	ss << date.month << "-";
	if(date.day < 10)
		ss << "0";
	ss << date.day;
	return os << ss.str();
}



std::ostream& operator<<(std::ostream& os,const DataSheet::Time& time)
{
	std::ostringstream ss; // setw should work
	ss << time.hour << ":";
	if (time.minute < 10)
		ss << "0";
	ss << time.minute;
	return os << ss.str();
}



template<typename T>
optional<T> parse_optional(T(*f)(const std::string&),const std::string& str)
{
	if (str.empty())
		return optional<T>();
	return optional<T>(f(str));
}




short parse_short(const std::string& str)
{
	int i = parse_int(str);
	if (i<std::numeric_limits<short>::min() || i>std::numeric_limits<short>::max())
		throw std::out_of_range("value out of range");
	return short(i);
}




DataSheet::Time parse_time(const std::string& timestr)
{
	auto tokens = split(timestr,':');
	if (tokens.size()!=2)
		throw std::invalid_argument("not a time");
	short h = parse_short(tokens.at(0));
	short m = parse_short(tokens.at(1));
	if (h<0 || h>23)
		throw std::out_of_range("hours out of range");
	if (m<0 || m>59)
		throw std::out_of_range("minutes out of range");
	return {h,m};
}



DataSheet::Date parse_date(const std::string& datestr)
{
	auto date_tokens = split(datestr,'-');
	if (date_tokens.size()!=3)
		throw std::invalid_argument("not a ISO date");

	int year = parse_int(date_tokens.at(0));	
	short month = parse_short(date_tokens.at(1));
	short day = parse_short(date_tokens.at(2));
	return {year,month,day};
}



std::string parse_text(const std::string& text)
{
	if (text.size()<2 || text.front()!='"' || text.back()!='"')
		throw std::invalid_argument("expected quoted text");

	return text.substr(1,text.size()-2);
}



bool parse_sign(const std::string& sign)
{
	if (sign=="x")
		return true;
	if (sign=="")
		return false;
	throw std::invalid_argument("expected 'x' or nothing");
}



bool parse_adjust(const std::string& str)
{
	return str.find('a')!=std::string::npos;
}

bool parse_plus(const std::string& str)
{
	return str.find('+')!=std::string::npos;
}

bool parse_minus(const std::string& str)
{
	return str.find('-')!=std::string::npos;
}

bool parse_reverse(const std::string& str)
{
	return str.find('r')!=std::string::npos;
}

bool parse_reverse_disabled(const std::string& str)
{
	return str.find('R')!=std::string::npos;
}

bool parse_invert(const std::string& str)
{
	return str.find('i')!=std::string::npos;
}

bool parse_bad(const std::string& str)
{
	return str.find('b')!=std::string::npos;
}

bool parse_cancel(const std::string& str)
{
	return str.find('c')!=std::string::npos;
}

bool parse_cancel_disabled(const std::string& str)
{
	return str.find('C')!=std::string::npos;
}

bool parse_inverted(const std::string& str)
{
	return str.find('i')!=std::string::npos;
}

bool parse_desk_in_sw(const std::string& str)
{
	return str.find('s')!=std::string::npos;
}

bool parse_desk_in_zw(const std::string& str)
{
	return str.find('z')!=std::string::npos;
}

DataSheet::WeatherObscurations
parse_weather_obscurations(const std::string& str)
{
	if (str.find('f')!=std::string::npos)
		return DataSheet::WeatherObscurations::Fog;
	if (str.find('m')!=std::string::npos)
		return DataSheet::WeatherObscurations::Mist;
	if (str.find('h')!=std::string::npos)
		return DataSheet::WeatherObscurations::Haze;

	return DataSheet::WeatherObscurations::Clear;
}

bool parse_weather_clouds(const std::string& str)
{
	return str.find('c')!=std::string::npos;
}

bool parse_weather_rain(const std::string& str)
{
	return str.find('r')!=std::string::npos;
}

bool parse_weather_wind(const std::string& str)
{
	return str.find('w')!=std::string::npos;
}

bool parse_desk_disabled(const std::string& str)
{
	return str.find('D')!=std::string::npos;
}

bool parse_visitors(const std::string& str)
{
	return str.find('v')!=std::string::npos;
}

bool parse_corrugated_paper(const std::string& str)
{
	return str.find('p')!=std::string::npos;
}



void parse_header_line(const std::string& line,DataSheet& ds)
{	
	int column = 0;
	try {
		auto tokens = split(line,CSV_DELIMITER);
		if (tokens.size()!=9) {
			throw ParseException{"expected 9 columns",0};
		}
		column++;
		ds.no = parse_int(tokens.at(0));
		column++;
		ds.date = parse_date(tokens.at(1));
		column++;		
		ds.mean_observation_time = parse_time(tokens.at(2));
		column++;
		ds.local_mean_time = parse_optional(parse_time,tokens.at(3));
		column++;
		ds.sidereal_time = parse_time(tokens.at(4));
		column++;
		ds.weight = parse_optional(parse_double,tokens.at(5));
		column++;
		ds.fringes = parse_optional(parse_int,tokens.at(6));
		column++;
		ds.sign_correct = parse_sign(tokens.at(7));
		column++;
		ds.comment = parse_text(tokens.at(8));


	}
	catch(std::exception& e) {
		throw ParseException {e.what(),column};
	}
}



optional<DataSheet::Thermometers>
parse_thermometer_line(const std::string& line,optional<DataSheet::Time>& thermometer_time)
{
	auto tokens = split(line,CSV_DELIMITER);
	if (tokens.size()<5)
		throw ParseException{"expected 5 columns",0};

	thermometer_time = parse_optional(parse_time,tokens.at(0));

	auto N = parse_optional(parse_double,tokens.at(1));
	auto E = parse_optional(parse_double,tokens.at(2));
	auto S = parse_optional(parse_double,tokens.at(3));
	auto W = parse_optional(parse_double,tokens.at(4));

	if (!N || !E || !S || !W)
		return {};


	return DataSheet::Thermometers {*N,*E,*S,*W};
}



void parse_weather(const std::string& line,DataSheet& ds)
{
	auto tokens = split(line,CSV_DELIMITER);
	if (tokens.size()<5)
		throw ParseException{"expected 5 or 6 columns",0};
	if (tokens.size()<6)
		return; // no weather codes

	ds.weather_obscuration = parse_weather_obscurations(tokens.at(5));
	ds.weather_clouds = parse_weather_clouds(tokens.at(5));
	ds.weather_rain = parse_weather_rain(tokens.at(5));
	if (ds.weather_rain)
		ds.weather_clouds = true;
	ds.weather_wind = parse_weather_wind(tokens.at(5));
}


void parse_start_end_times(const std::string& line,const Options& options,DataSheet& ds)
{
	auto tokens = split(line,CSV_DELIMITER);
	if (tokens.size()<2 || tokens.size()>3)
		throw ParseException {"expected 2 or 3 columns",0};

	ds.start_time = parse_time(tokens.at(0));
	ds.end_time = parse_time(tokens.at(1));
	if (tokens.size() == 3) {
		ds.inverted = !options.ignore_inverted_theory && parse_inverted(tokens.at(2));
		ds.reverse = !options.ignore_reverse_sheet && parse_reverse(tokens.at(2));
		ds.desk_in_sw = parse_desk_in_sw(tokens.at(2));
		if (!options.ignore_sw)
			ds.desk_in_sw = ds.desk_in_sw || parse_desk_in_zw(tokens.at(2));
		ds.desk_disabled = parse_desk_disabled(tokens.at(2));
		ds.visitors = parse_visitors(tokens.at(2));
		ds.corrugated_paper = parse_corrugated_paper(tokens.at(2));
	}
}



void parse_data_line(const std::string& line,const Options& options,DataSheet& ds)
{
	auto tokens = split(line,CSV_DELIMITER);
	if (tokens.size()<17 || tokens.size()>18)
		throw ParseException{"expected 17 or 18 columns",0};

	DataSheet::Turn turn;
	if (tokens.size()==18) {
		turn.adjust = parse_adjust(tokens.at(17));
		turn.bad = !options.ignore_bad && parse_bad(tokens.at(17));
		if (!options.ignore_cancel) {
			turn.cancel = parse_cancel(tokens.at(17));
			if (options.ignore_cancel_disabled)
				turn.cancel = turn.cancel || parse_cancel_disabled(tokens.at(17));
		}
		if (!options.ignore_reverse) {
			turn.reverse = parse_reverse(tokens.at(17));
			if (options.ignore_reverse_disabled)
				turn.reverse = turn.reverse || parse_reverse_disabled(tokens.at(17));
		}
		turn.invert = !options.ignore_invert && parse_invert(tokens.at(17));
		turn.plus = parse_plus(tokens.at(17));
		turn.minus = parse_minus(tokens.at(17));
		tokens.pop_back();		
	}

	int i=0;
	try {
		for (const std::string& token : tokens) {
			turn.distances.at(i) = parse_short(token);
			i++;
		}
	}
	catch(std::exception& e) {
		throw ParseException{e.what(),i+1};
	}

	ds.turns.push_back(std::move(turn));
}



void parse_line(const std::string& line,int line_no,const Options& options,DataSheet& ds)
{
	if (line_no <= HEADER_LINES) { // header
		switch(line_no) {
		case 1:
			parse_header_line(line,ds);
			break;
		case 2:
			ds.thermometers_start = parse_thermometer_line(line,ds.thermometers_start_time);
			parse_weather(line,ds);
			break;
		case 3:
			ds.thermometers_end = parse_thermometer_line(line,ds.thermometers_end_time);
			break;
		case 4:
			parse_start_end_times(line,options,ds);
			break;
		default:
			throw std::runtime_error(SSSTR("unknown line no " << line_no));
		}
	}
	else { // data
		parse_data_line(line,options,ds);
	}

}



bool validate_filename(const std::string& filename, const DataSheet& data_sheet)
{
	std::stringstream name;
	name << "dcm_" << data_sheet.date.year << "-";
	if (data_sheet.date.month<10)
		name << "0";
	name << data_sheet.date.month << "-";
	if (data_sheet.date.day<10)
		name << "0";
	name <<  data_sheet.date.day;
	name << "_" << data_sheet.no << ".csv";
	return ends_with(filename,name.str());
}



DataSheet load_data_sheet_csv(const std::string& filename,const Options& options)
{
	DataSheet ds;


	std::ifstream fs(filename);
	if (!fs) {
		std::cerr << "Error opening file " << filename << "\n";
		throw ExitException();
	}

	std::string line;
	int line_no = 0;
	while (!fs.eof()) {
		std::getline(fs,line);
		if (fs.fail() && !fs.eof()) {
			std::cerr << "Error reading file " << filename << "\n";
			throw ExitException();
		}

		if (!line.empty()) {
			line_no++;
			try {
				parse_line(line,line_no,options,ds);
			}
			catch(ParseException& e) {
				std::cerr << "Error at line " << line_no;
				if (e.pos>0)
					std::cerr << " column " << e.pos;
				std::cerr << " while parsing file " << filename << ": ";
				std::cerr << e.message << "\n";
				throw ExitException();
			}
			catch(std::exception& e) {
				std::cerr << "Error in line " << line_no << " while parsing file " << filename << ": ";
				std::cerr << e.what() << "\n";
				throw ExitException();
			}
		}

	}

	if (line_no <= HEADER_LINES) {
		std::cerr << "Invalid file format in file " << filename << "\n";
		throw ExitException();
	}	


	//if (!validate_filename(filename,ds)) {
		//std::cerr << "warning: filename does not match for file " << filename << "\n";
	//}

	return ds;
}




std::string data_sheet_csv_sign_value(bool sign)
{
	return sign ? "x" : "";
}






std::string data_sheet_attributes(const DataSheet& data_sheet)
{
	std::string attributes;
	if (data_sheet.desk_disabled)
		attributes += "D";
	if (data_sheet.desk_in_sw)
		attributes += "s";
	if (data_sheet.inverted)
		attributes += "i";
	if (data_sheet.reverse)
		attributes += "r";
	if (data_sheet.visitors)
		attributes += "v";
	if (data_sheet.corrugated_paper)
		attributes += "p";

	return attributes;
}


std::string data_sheet_weather(const DataSheet& data_sheet)
{
	std::string codes;
	switch(data_sheet.weather_obscuration)
	{
	case DataSheet::WeatherObscurations::Clear:
		break;
	case DataSheet::WeatherObscurations::Haze:
		codes += "h";
		break;
	case DataSheet::WeatherObscurations::Mist:
		codes += "m";
		break;
	case DataSheet::WeatherObscurations::Fog:
		codes += "f";
		break;
	default:
		throw std::runtime_error("unknown atmospheric obscuration");
	}

	if (data_sheet.weather_clouds)
		codes += "c";
	if (data_sheet.weather_rain)
		codes += "r";
	if (data_sheet.weather_wind)
		codes += "w";


	return codes;
}





void write_data_sheet_csv(std::ostream& os,const DataSheet& data_sheet)
{

	output_args_separated(os,CSV_DELIMITER,data_sheet.no,data_sheet.date,
					 data_sheet.mean_observation_time,data_sheet.local_mean_time,data_sheet.sidereal_time,
					 data_sheet.weight,data_sheet.fringes,data_sheet_csv_sign_value(data_sheet.sign_correct),
					 quote(data_sheet.comment));
	os << "\n";

	output_args_separated(os,CSV_DELIMITER,data_sheet.thermometers_start_time,"");
	if (data_sheet.thermometers_start.has_value()) {
		output_args_separated(os,CSV_DELIMITER,
							  data_sheet.thermometers_start->N,
							  data_sheet.thermometers_start->E,
							  data_sheet.thermometers_start->S,
							  data_sheet.thermometers_start->W,
							  data_sheet_weather(data_sheet));
	}
	os << "\n";

	output_args_separated(os,CSV_DELIMITER,data_sheet.thermometers_end_time,"");
	if (data_sheet.thermometers_end.has_value()) {
		output_args_separated(os,CSV_DELIMITER,
							  data_sheet.thermometers_end->N,
							  data_sheet.thermometers_end->E,
							  data_sheet.thermometers_end->S,
							  data_sheet.thermometers_end->W);
	}
	os << "\n";


	output_args_separated(os,CSV_DELIMITER,data_sheet.start_time,data_sheet.end_time,
						  data_sheet_attributes(data_sheet));
	os << "\n";

	for (const auto& turn : data_sheet.turns) {
		output_separated(os,turn.distances,CSV_DELIMITER); // TODO to short?
		os << "\n";
	}

	os << "\n";
}



double time_to_h(const DataSheet::Time& t)
{
	return t.hour + t.minute/60.0;
}



DataSheet::Time h_to_time(double h)
{
	h = std::fmod(h,24);
	double hh=0;
	auto hm = std::modf(h,&hh);
	return {short(hh),short(std::min(59.0,std::ceil(60*hm)))}; // bad, but good enough
}



double max_T(const DataSheet::Thermometers& th)
{
	double max = th.N;
	if (th.E>max)
		max = th.E;
	if (th.S>max)
		max = th.S;
	if (th.W>max)
		max = th.W;

	return max;
}



double min_T(const DataSheet::Thermometers& th)
{
	double min = th.N;
	if (th.E<min)
		min = th.E;
	if (th.S<min)
		min = th.S;
	if (th.W<min)
		min = th.W;

	return min;
}




bool has_near_duplicates(const DataSheet& data_sheet,int turns,size_t& turn1, size_t& turn2)
{
	for (size_t i=0;i<data_sheet.turns.size();i++) {
		for (size_t k=i+1;k<=i+turns && k<data_sheet.turns.size();k++) {
			if (std::equal(data_sheet.turns.at(i).distances.begin(),std::prev(data_sheet.turns.at(i).distances.end()),
					   data_sheet.turns.at(k).distances.begin())) {
				turn1 = i+1;
				turn2 = k+1;
				return true;
			}
		}
	}

	return false;
}



int minutes(const DataSheet::Time& t)
{
	return t.hour*60 + t.minute;
}



bool validate_start_end_times(const DataSheet& data_sheet,std::ostream& os)
{
	int duration;
	if (data_sheet.start_time.hour <= data_sheet.end_time.hour) {
		duration = minutes(data_sheet.end_time)-minutes(data_sheet.start_time);
	}
	else {
		if (data_sheet.end_time.hour!=0) {
			os << "error: start time > end time\n";
			return false;
		}
		duration = minutes(data_sheet.end_time) + 24*60-minutes(data_sheet.start_time);
	}

	if (duration<0) {
		os << "error: start time > end time\n";
		return false;
	}

	auto mind = MIN_DURATION/20.0*data_sheet.turns.size();
	auto maxd = MAX_DURATION/20.0*data_sheet.turns.size();
	if (duration<mind) {
		os << "warning: start-end duration less than " << mind << " min\n";
		return false;
	}
	else if (duration>maxd) {
		os << "warning: start-end duration greater than " << maxd << " min\n";
		return false;
	}

	return true;
}



int duration_in_min(const DataSheet::Time& t1,const DataSheet::Time& t2)
{
	int duration;

	duration = minutes(t2)-minutes(t1);
	if (t1.hour > t2.hour)
		duration += 24*60;

	return duration;
}



bool validate_mean_times(const DataSheet& data_sheet,std::ostream& os)
{
	bool ok = true;
	double mean = duration_in_min(data_sheet.start_time,data_sheet.end_time)/2.0;
	double mean_duration = duration_in_min(data_sheet.start_time,data_sheet.mean_observation_time);
	if (!approximate_delta(mean,mean_duration,2)) {
		os << "error: mean observation time does not match\n";
		ok = false;
	}

	return ok;
}



Calendar calendar_date(const DataSheet& data_sheet,bool UT)
{
	Calendar calendar {data_sheet.date.year, data_sheet.date.month, double(data_sheet.date.day)};
	if (UT && data_sheet.local_mean_time.has_value()) {
		calendar.day += time_to_h(data_sheet.local_mean_time.value())/24.0;
	}
	else {
		calendar.day += time_to_h(data_sheet.mean_observation_time)/24.0;
	}
	return calendar;
}



bool validate_sidereal_time(const DataSheet& data_sheet,std::ostream& os)
{
	// Millers sidereal time is constantly ~3 min before the calculated value
	const double hprec = 1.0/60.0*3;

	Calendar calendar = calendar_date(data_sheet,true);
	auto theta = sidereal_time(calendar,0);
	if (!approximate_delta_periodic(rad_to_h(theta),time_to_h(data_sheet.sidereal_time),hprec,24)) {
		os << "error: sidereal time " << data_sheet.sidereal_time;
		os << " does not match calculated value " << h_to_time(rad_to_h(theta)) << "\n";
		return false;
	}

	return true;
}






bool validate_temperatur_difference(const optional<DataSheet::Thermometers>& thermometers,std::ostream& os)
{	
	if (thermometers.has_value()) {
		auto minT = min_T(*thermometers);
		auto maxT = max_T(*thermometers);
		auto delta = std::abs(maxT-minT); // degree Celsius
		if (delta > 2) {
			os << "warning: temperature difference " << delta << "\n";
			return false;
		}
	}

	return true;
}



bool validate_temperature_change(double t1,double t2,std::ostream& os)
{
	auto dt = std::abs(t2-t1);
	if (dt > 2) {
		os << "warning: temperature change of " << dt << "\n";
		return false;
	}
	return true;
}



bool validate_temperature_change(const optional<DataSheet::Thermometers>& t1,
								 const optional<DataSheet::Thermometers>& t2,std::ostream& os)
{
	if (!t1 || !t2)
		return true;

	bool ok = true;
	ok = ok && validate_temperature_change(t1->N,t2->N,os);
	ok = ok && validate_temperature_change(t1->E,t2->E,os);
	ok = ok && validate_temperature_change(t1->S,t2->S,os);
	ok = ok && validate_temperature_change(t1->W,t2->W,os);

	return ok;
}



bool validate_temperatures(const DataSheet& data_sheet,std::ostream& os)
{
	bool ok = true;
	ok = ok && validate_temperatur_difference(data_sheet.thermometers_start,os);
	ok = ok && validate_temperatur_difference(data_sheet.thermometers_end,os);
	ok = ok && validate_temperature_change(data_sheet.thermometers_start,data_sheet.thermometers_end,os);
	return ok;
}



bool validate(const DataSheet& data_sheet,const std::string& filename,std::ostream& os)
{
	bool ok = true;

	ok = ok && validate_start_end_times(data_sheet,os);
	ok = ok && validate_mean_times(data_sheet,os);
	ok = ok && validate_sidereal_time(data_sheet,os);
	ok = ok && validate_temperatures(data_sheet,os);

	if (data_sheet.turns.size() < USUAL_NUMBER_OF_TURNS) {
		os << "warning: less than " << USUAL_NUMBER_OF_TURNS << " data rows\n";
		ok = false;
	}

	size_t turn1,turn2;
	if (has_near_duplicates(data_sheet,4,turn1,turn2)) {
		os << "warning: found duplicate turns " << turn1 << " and " << turn2 << "\n";
		ok = false;
	}

	if (!ok)
		os << "  in file " << filename << "\n";

	return ok;
}



double mean_T(const DataSheet::Thermometers& th)
{
	return (th.N + th.E + th.S + th.W)/4.0;
}



optional<double> mean_T(const DataSheet& data_sheet)
{	
	if (!(data_sheet.thermometers_start || data_sheet.thermometers_end))
		return {};
	
	double T = 0.0;	
	int n=0;
	if (data_sheet.thermometers_start) {
		T += mean_T(*data_sheet.thermometers_start);
		n++;
	}
	if (data_sheet.thermometers_end) {
		T += mean_T(*data_sheet.thermometers_end);
		n++;
	}

	return T/n;
}




optional<double> max_TD(const DataSheet& data_sheet)
{
	if (!(data_sheet.thermometers_start || data_sheet.thermometers_end))
		return {};

	double TD = -std::numeric_limits<double>::infinity();
	if (data_sheet.thermometers_start) {
		TD = max_T(*data_sheet.thermometers_start)-min_T(*data_sheet.thermometers_start);
	}

	if (data_sheet.thermometers_end) {
		auto d = max_T(*data_sheet.thermometers_end)-min_T(*data_sheet.thermometers_end);
		if (d>TD)
			TD = d;
	}

	return TD;
}



optional<double> mean_TD(const optional<DataSheet::Thermometers>& th1,
						 const optional<DataSheet::Thermometers>& th2)
{
	if (!th1 && !th2)
		return {};

	int n=0;
	double TD = 0;
	if (th1) {
		TD += max_T(*th1)-min_T(*th1);
		n++;
	}

	if (th2) {
		TD += max_T(*th2)-min_T(*th2);
		n++;
	}

	return TD/n;
}




double mean_dT(const DataSheet::Thermometers& th1,const DataSheet::Thermometers& th2)
{
	return ((th2.N + th2.E + th2.S + th2.W) - (th1.N + th1.E + th1.S + th1.W))/4.0;
}



bool has_dT(const DataSheet& data_sheet)
{
	return data_sheet.thermometers_start.has_value() && data_sheet.thermometers_end.has_value();
}




optional<double> mean_dT(const DataSheet& data_sheet)
{	
	if (!has_dT(data_sheet)) 
		return {};
	
	auto dT = mean_dT(*data_sheet.thermometers_start,*data_sheet.thermometers_end);		
	
	DataSheet::Time t1,t2;
	if (data_sheet.thermometers_start_time) {
		t1 = *data_sheet.thermometers_start_time;			
	}
	else {
		t1 = data_sheet.start_time;			
	}
	if (data_sheet.thermometers_end_time) {
		t2 = *data_sheet.thermometers_end_time;				
	}
	else {
		t2 = data_sheet.end_time;				
	}
	
	auto dt = periodic_distance(time_to_h(t1),time_to_h(t2),24.); // handle midnight
	return dT/dt * 0.25; // normalise temp change to 15 mins	
}



double mean_dT(const DataSheet& prev,const DataSheet& next)
{
	double T1,T2;
					
	if (prev.thermometers_end.has_value()) {
		T1 = mean_T(*prev.thermometers_end);
	}
	else if (prev.thermometers_start.has_value()) {
		T1 = mean_T(*prev.thermometers_start);
	}
	else 
		return NAN;
	
	if (next.thermometers_start.has_value()) {
		T2 = mean_T(*next.thermometers_start);		
	}
	else if (next.thermometers_end.has_value()) {
		T2 = mean_T(*next.thermometers_end);		
	}
	else 
		return NAN;
	
	auto jd1 = julian_date(calendar_date(prev));	
	auto jd2 = julian_date(calendar_date(next));
	auto dT = T2-T1;
	auto dt = (jd2-jd1)*24;
	return dT/dt * 0.25; // normalise temp change to 15 mins
}




DataSheetStats data_sheet_stats(const DataSheet& data_sheet)
{
	DataSheetStats stats;

	double pos_drift = 0;
	double neg_drift = 0;
	int np=0;
	int nn=0;
	int plusminus = 0;
	for (const auto& turn : data_sheet.turns) {
		if (turn.adjust)
			stats.number_of_adjust++;

		if (!(turn.bad || turn.cancel))
			plusminus++;
		if (turn.plus)
			stats.number_of_plus += plusminus;
		if (turn.minus)
			stats.number_of_minus += plusminus;
		if (turn.plus || turn.minus)
			plusminus = 0;

		if (turn.bad || turn.cancel)
			continue;


		short drift = turn.distances.at(16)-turn.distances.at(0);
		//if (turn.invert || turn.reverse)
			//drift = -drift; dont transform?

		if (drift<0) {
			neg_drift += drift;
			nn++;
		}
		else {
			pos_drift += drift;
			np++;
		}

	}
	//stats.drift = (pos_drift+neg_drift)/(np+nn);
	stats.drift = (pos_drift*np + neg_drift*nn)/sqr(np+nn);
	stats.abs_drift = (pos_drift - neg_drift)/(np+nn);
	stats.mean_T = mean_T(data_sheet);
	stats.max_TD = max_TD(data_sheet);	
	stats.mean_dT = mean_dT(data_sheet);

	return stats;
}




Epoch epoch(short month)
{
	switch(month) {
	case 2:
		return Epoch::Feb;
	case 3:
	case 4:
		return Epoch::Apr;
	case 7:
	case 8:
		return Epoch::Aug;
	case 9:
		return Epoch::Sep;
	default:
		throw std::runtime_error("unsupported month for epoch");
	}
}

Epoch epoch(const DataSheet& data_sheet)
{
	return epoch(data_sheet.date.month);
}


std::ostream& operator << (std::ostream& os,const Epoch& epoch)
{
	switch(epoch) {
	case Epoch::Feb:
		os << "feb";		
		break;
	case Epoch::Apr:		
		os << "apr";		
		break;
	case Epoch::Aug:
		os << "aug";		
		break;
	case Epoch::Sep:
		os << "sep";		
		break;		
	default:
		os << "???";
	}
	return os;
}


TimeOfDay time_of_day(const DataSheet& data_sheet)
{
	auto time = time_to_h(data_sheet.mean_observation_time);

	// values from calsky.com
	switch(epoch(data_sheet)) {
	case Epoch::Feb:
		// reference date 1926-02-07, according to Miller (sheet 91) "sunrise about 6:40"
		return time_of_day(time_to_h({6,45}),time_to_h({17,28}),time);		
	case Epoch::Apr:
		// reference date 1925-04-03
		return time_of_day(time_to_h({5,37}),time_to_h({18,14}),time);		
	case Epoch::Aug:
		// reference date 1925-07-30, according to Miller (sheet 26) "sunrise about 5:10am"
		return time_of_day(time_to_h({5,02}),time_to_h({18,55}),time);		
	case Epoch::Sep:
		// reference date 1925-09-16
		return time_of_day(time_to_h({5,36}),time_to_h({17,58}),time);		
	default:
		throw std::runtime_error("unsupported epoch for time of day");
	}
}




IntegerInterval month_interval(const DataSheet& data_sheet)
{
	auto ep = epoch(data_sheet);
	switch(ep) {
	case Epoch::Apr:
		return {3,4};
	case Epoch::Aug:
		return {7,8};
	case Epoch::Sep:
		return {9,9};
	case Epoch::Feb:
		return {2,2};
	default:
		return {0,0};
	}
}



}//aether
