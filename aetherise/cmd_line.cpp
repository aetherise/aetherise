#include "cmd_line.h"

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>

namespace aether {


bool equal(const char* a,const char* b)
{
	if (a==nullptr || b==nullptr)
		return false;
	return std::strcmp(a,b) == 0;
}



Filter::Interval parse_interval(int argc,char* argv[],int& i)
{
	Filter::Interval range;

	i++;
	if (i<argc) {
		std::string rstr = argv[i];
		if (rstr.size()<3 || rstr.front()!='[' || rstr.back()!=']') {
			std::cerr << "Invalid interval format for option " << argv[i-1] << "\n";
			throw ExitException();
		}

		rstr = rstr.substr(1,rstr.size()-2);
		auto tokens = split(rstr,',');
		if (tokens.size()!=2) {
			std::cerr << "Invalid interval " << argv[i] << "\n";
			throw ExitException();
		}

		try {
			if (!tokens.at(0).empty())
				range.min = std::stod(tokens.at(0));
			if (!tokens.at(1).empty())
				range.max = std::stod(tokens.at(1));
		}
		catch(std::exception& ) {
			std::cerr << "Invalid value in interval " << argv[i] << "\n";
			throw ExitException();
		}

		if (range.max < range.min) {
			std::cerr << "Empty interval " << argv[i] << "\n";
			throw ExitException();
		}

	}
	else {
		std::cerr << "Missing interval for option " << argv[i-1] << "\n";
		throw ExitException();
	}

	return range;
}


void parse_reduction_argument(int argc,char* argv[],int& i,Options& options)
{

	i++;
	if (i<argc) {
		std::string argstr = argv[i];
		if (argstr == "Miller") {
			options.reduction_method = Options::DataReductionMethod::Miller;
		}
		else if (argstr == "separate") {
			options.reduction_method = Options::DataReductionMethod::Separate;
		}
		else {
			std::cerr << "Unknown reduction method " << argstr << "\n";
			throw ExitException();
		}
	}
	else {
		std::cerr << "Missing argument for option " << argv[i-1] << "\n";
		throw ExitException();
	}
}

void parse_theory_argument(int argc,char* argv[],int& i,Options& options)
{

	i++;
	if (i<argc) {
		std::string argstr = argv[i];
		if (argstr == "classic") {
			options.theory = Options::Theory::Classic;
		}	
		else if (argstr == "aether") {
			options.theory = Options::Theory::Aether;
		}
		else if (argstr == "relativity") {
			options.theory = Options::Theory::Relativity;
		}
		else {
			std::cerr << "Unknown theory " << argstr << "\n";
			throw ExitException();
		}
	}
	else {
		std::cerr << "Missing argument for option " << argv[i-1] << "\n";
		throw ExitException();
	}

}


void parse_minimizer_argument(int argc,char* argv[],int& i,Options& options)
{

	i++;
	if (i<argc) {
		std::string argstr = argv[i];
		if (argstr == "grad") {
			options.minimizer = Options::Minimizer::Grad;
		}
#ifdef AETHER_MINUIT
		else if (argstr == "Minuit2") {
			options.minimizer = Options::Minimizer::Minuit2;
		}
#endif
		else {
			std::cerr << "Unknown minimization method " << argstr << "\n";
			throw ExitException();
		}
	}
	else {
		std::cerr << "Missing argument for option " << argv[i-1] << "\n";
		throw ExitException();
	}

}

void parse_aggregate_argument(int argc,char* argv[],int& i,Options& options)
{

	i++;
	if (i<argc) {
		std::string argstr = argv[i];
		if (argstr == "mean") {
			options.aggregation_method = Options::AggregationMethod::Mean;
		}
		else if (argstr == "sidereal") {
			options.aggregation_method = Options::AggregationMethod::Sidereal;
		}
		else if (argstr == "diff_chi2") {
			options.aggregation_method = Options::AggregationMethod::DiffChi;
		}
		else if (argstr == "params") {
			options.aggregation_method = Options::AggregationMethod::Params;
		}
		else if (argstr == "model_chi") {
			options.aggregation_method = Options::AggregationMethod::ModelChi;
		}
		else if (argstr == "fit") {
			options.aggregation_method = Options::AggregationMethod::Fit;
		}
		else if (argstr == "list") {
			options.aggregation_method = Options::AggregationMethod::List;
		}
		else if (argstr == "test") {
			options.aggregation_method = Options::AggregationMethod::Test;
		}
		else if (argstr == "signals") {
			options.aggregation_method = Options::AggregationMethod::Signals;
		}
		else {
			std::cerr << "Unknown aggregation method " << argstr << "\n";
			throw ExitException();
		}
	}
	else {
		std::cerr << "Missing argument for option " << argv[i-1] << "\n";
		throw ExitException();
	}

}


void parse_data_file_line(const std::string& line,std::vector<double>& values)
{
	auto str_values = split(line,';');
	if (str_values.size()<16) {
		std::cerr << "Data file line must contain at least 16 values\n";
		throw ExitException();
	}

	for (auto& str : str_values) {
		values.push_back(parse_double(str)); // throws
	}
}


void read_data_file_line(std::ifstream& fs,const std::string& filename,
						 std::string& line,std::vector<double>& values)
{
	if (!fs.eof()) {
		std::getline(fs,line);
		if (!fs.eof()) { // double check eof to allow empty line
			if (fs.fail()) {
				std::cerr << "Error reading file " << filename << "\n";
				throw ExitException();
			}
			parse_data_file_line(line,values);
		}
	}
}


void load_data_file(const std::string& filename,Options& options)
{
	std::ifstream fs(filename);
	if (!fs) {
		std::cerr << "Error opening file " << filename << "\n";
		throw ExitException();
	}

	std::string line;
	try {
		read_data_file_line(fs,filename,line,options.data1);
		read_data_file_line(fs,filename,line,options.data2);
		read_data_file_line(fs,filename,line,options.data3);
		read_data_file_line(fs,filename,line,options.data4);
	}
	catch(std::invalid_argument& e) {
		std::cerr << "Error parsing file " << filename << ": " << e.what() << "\n";
		throw ExitException();
	}
}



void parse_data_filename(int argc,char* argv[],int& i,Options& options)
{
	i++;
	if (i<argc) {
		options.data_filename = argv[i];
		load_data_file(options.data_filename,options);
	}
	else {
		std::cerr << "Missing argument for option " << argv[i-1] << "\n";
		throw ExitException();
	}
}



bool are_valid_ignore_codes(const std::string& codes)
{
	if (codes=="all" || codes=="all!")
		return true;

	const char* valid_codes = "-ibrcRCz";
	for (const char c : codes ) {
		auto p = std::strchr(valid_codes,c);
		if (p == nullptr)
			return false;
	}
	return true;
}



void parse_ignore_codes(char* argv,Options& options)
{
	std::string codes(argv);

	if (!are_valid_ignore_codes(codes)) {
		std::cerr << "ERROR: invalid ignore code\n";
		throw ExitException();
	}


	if (codes == "all" || codes == "all!") {
		if (codes=="all!") {
			options.ignore_reverse = true;
			options.ignore_cancel = true;
		}
		options.ignore_reverse_sheet = true;
		options.ignore_invert = true;
		options.ignore_bad = true;
		options.ignore_inverted_theory = true;
		options.ignore_reverse_disabled = true;
		options.ignore_cancel_disabled = true;
		options.ignore_sw = true;
	}	
	else {
		if (codes.find('-')!=std::string::npos)
			options.ignore_reverse_sheet = true;
		if (codes.find('r')!=std::string::npos)
			options.ignore_reverse = true;
		if (codes.find('c')!=std::string::npos)
			options.ignore_cancel = true;
		if (codes.find('i')!=std::string::npos)
			options.ignore_invert = true;
		if (codes.find('b')!=std::string::npos)
			options.ignore_bad = true;
		if (codes.find('R')!=std::string::npos)
			options.ignore_reverse_disabled = true;
		if (codes.find('C')!=std::string::npos)
			options.ignore_cancel_disabled = true;
		if (codes.find('z')!=std::string::npos)
			options.ignore_sw = true;

	}
}


void parse_ignore(int argc,char* argv[],int& i,Options& options)
{
	i++;
	if (i<argc) {
		parse_ignore_codes(argv[i],options);
	}
	else {
		std::cerr << "Missing argument for option " << argv[i-1] << "\n";
		throw ExitException();
	}
}



template<typename F>
void parse_numeric_argument(int argc,char* argv[],int& i,F f)
{
	i++;
	if (i<argc) {
		try {
			f();
		} catch(std::invalid_argument& e) {
			std::cerr << "Invalid value for option " << argv[i-1] << ": " << e.what() << "\n";
			throw ExitException();
		}
	}
	else {
		std::cerr << "Missing argument for option " << argv[i-1] << "\n";
		throw ExitException();
	}
}


void parse_signals_dTD_argument(int argc,char* argv[],int& i,Options& options)
{
	parse_numeric_argument(argc,argv,i,[&](){
		options.signals_dTD = parse_double(argv[i]);
		if (options.signals_dTD<=0)	{
			std::cerr << "Expected value >0 for option " << argv[i-1] << "\n";
			throw ExitException();
		}	
	});
}



void parse_signals_ddT_argument(int argc,char* argv[],int& i,Options& options)
{
	parse_numeric_argument(argc,argv,i,[&](){
		options.signals_ddT = parse_double(argv[i]);		
		if (options.signals_ddT<=0)	{
			std::cerr << "Expected value >0 for option " << argv[i-1] << "\n";
			throw ExitException();
		}	
	});
}




void parse_signals_dt_argument(int argc,char* argv[],int& i,Options& options)
{
	parse_numeric_argument(argc,argv,i,[&](){
		options.signals_dt = parse_double(argv[i]);		
		if (options.signals_dt<=0)	{
			std::cerr << "Expected value >0 for option " << argv[i-1] << "\n";
			throw ExitException();
		}	
	});
}



void parse_sim_seed_argument(int argc,char* argv[],int& i,Options& options)
{
	parse_numeric_argument(argc,argv,i,[&](){
		auto ul = parse_ulong(argv[i]);
		if (ul > std::numeric_limits<unsigned int>::max())	{
			std::cerr << "Value out of range for option " << argv[i-1] << "\n";
			throw ExitException();
		}	
		options.sim_seed = unsigned(ul);		
	});
}



void parse_delta_chi_squared_argument(int argc,char* argv[],int& i,Options& options)
{
	i++;
	if (i<argc) {
		try {
			options.delta_chi_squared = parse_double(argv[i]);			
			if (options.delta_chi_squared<1) {
				std::cerr << "Expected value >=1 for option " << argv[i-1] << "\n";
				throw ExitException();
			}
		} catch(std::invalid_argument& e) {
			std::cerr << "Invalid value for option " << argv[i-1] << ": " << e.what() << "\n";
			throw ExitException();
		}
	}
	else {
		std::cerr << "Missing argument for option " << argv[i-1] << "\n";
		throw ExitException();
	}
}


void parse_chi_squared_scale_argument(int argc,char* argv[],int& i,Options& options)
{
	i++;
	if (i<argc) {
		try {
			options.chi_squared_scale = parse_double(argv[i]);
			if (options.chi_squared_scale<=0) {
				std::cerr << "Expected positive value for option " << argv[i-1] << "\n";
				throw ExitException();
			}
		} catch(std::invalid_argument& e) {
			std::cerr << "Invalid value for option " << argv[i-1] << ": " << e.what() << "\n";
			throw ExitException();
		}
	}
	else {
		std::cerr << "Missing argument for option " << argv[i-1] << "\n";
		throw ExitException();
	}
}


TheoryParameters parse_params_argument(const char* arg)
{
	auto token = split(std::string(arg),',');
	if (token.size()!=3) {
		throw std::invalid_argument("expected 3 parameters");
	}

	TheoryParameters params;
	params.v = parse_double(token.at(0));
	params.a = h_to_rad(parse_double(token.at(1)));
	params.d = rad(parse_double(token.at(2)));
	return params;
}


void parse_theory_params_argument(int argc,char* argv[],int& i,Options& options)
{
	i++;
	if (i<argc) {
		try {
			options.theory_params = parse_params_argument(argv[i]);
		}
		catch(std::invalid_argument& e) {
			std::cerr << "Invalid argument for option " << argv[i-1] << ": " << e.what() << "\n";
			throw ExitException();
		}
	}
	else {
		std::cerr << "Missing argument for option " << argv[i-1] << "\n";
		throw ExitException();
	}
}


void parse_start_params_argument(int argc,char* argv[],int& i,Options& options)
{
	i++;
	if (i<argc) {
		try {
			options.start_params = parse_params_argument(argv[i]);
		}
		catch(std::invalid_argument& e) {
			std::cerr << "Invalid argument for option " << argv[i-1] << ": " << e.what() << "\n";
			throw ExitException();
		}
	}
	else {
		std::cerr << "Missing argument for option " << argv[i-1] << "\n";
		throw ExitException();
	}
}



void parse_index_of_refraction_argument(int argc,char* argv[],int& i,Options& options)
{
	i++;
	if (i<argc) {
		try {
			options.index_of_refraction = parse_double(argv[i]);

			if (options.index_of_refraction<1 || options.index_of_refraction>2) {
				std::cerr << "ERROR: index of refraction not in the valid interval [1, 2]\n";
				throw ExitException();
			}

			if (options.index_of_refraction<1.0002 || options.index_of_refraction>1.0003) {
				std::cerr << "WARNING: index of refraction not in the expected interval [1.0002, 1.0003]\n";
			}

		}
		catch(std::invalid_argument& e) {
			std::cerr << "Invalid argument for option " << argv[i-1] << ": " << e.what() << "\n";
			throw ExitException();
		}
	}
	else {
		std::cerr << "Missing argument for option " << argv[i-1] << "\n";
		throw ExitException();
	}
}



std::unordered_set<int> parse_fit_disable_argument(const char* argv)
{
	std::unordered_set<int> disabled_signals;
	
	std::string argvs(argv);
	auto nos = split(argvs,',');
	for (auto& no : nos) {
		disabled_signals.insert(parse_int(no));
	}
	
	return disabled_signals;
}



void parse_fit_disable(int argc,char* argv[],int& i,Options& options)
{
	i++;
	if (i<argc) {
		try {
			options.disabled_signals = parse_fit_disable_argument(argv[i]);
		}
		catch(std::invalid_argument& e) {
			std::cerr << "Invalid argument for option " << argv[i-1] << ": " << e.what() << "\n";
			throw ExitException();
		}
	}
	else {
		std::cerr << "Missing argument for option " << argv[i-1] << "\n";
		throw ExitException();
	}
}



void parse_option(int argc,char* argv[],int& i,Filter& filter,Action& action, bool& aggregate,Options& options)
{
	if (equal(argv[i],"-T")) {
		filter.T.push_back(parse_interval(argc,argv,i));
	}
	else if (equal(argv[i],"-dT")) {
		filter.dT.push_back(parse_interval(argc,argv,i));
	}
	else if (equal(argv[i],"-mean_dT")) {
		filter.mean_dT.push_back(parse_interval(argc,argv,i));
	}
	else if (equal(argv[i],"-TD")) {
		filter.TD.push_back(parse_interval(argc,argv,i));
	}
	else if (equal(argv[i],"-adjust")) {
		filter.adjust.push_back(parse_interval(argc,argv,i));
	}
	else if (equal(argv[i],"-no")) {
		filter.no.push_back(parse_interval(argc,argv,i));
	}
	else if (equal(argv[i],"-year")) {
		filter.year.push_back(parse_interval(argc,argv,i));
	}
	else if (equal(argv[i],"-month")) {
		filter.month.push_back(parse_interval(argc,argv,i));
	}
	else if (equal(argv[i],"-day")) {
		filter.day.push_back(parse_interval(argc,argv,i));
	}
	else if (equal(argv[i],"-time")) {
		filter.time.push_back(parse_interval(argc,argv,i));
	}
	else if (equal(argv[i],"-sidereal")) {
		filter.sidereal_time.push_back(parse_interval(argc,argv,i));
	}
	else if (equal(argv[i],"-weight")) {
		filter.weight.push_back(parse_interval(argc,argv,i));
	}
	else if (equal(argv[i],"-fringes")) {
		filter.fringes.push_back(parse_interval(argc,argv,i));
	}
	else if (equal(argv[i],"-amplitude")) {
		filter.amplitude.push_back(parse_interval(argc,argv,i));
	}
	else if (equal(argv[i],"-uncertainty")) {
		filter.uncertainty.push_back(parse_interval(argc,argv,i));
	}
	else if (equal(argv[i],"-theory_amp")) {
		filter.theory_amplitude.push_back(parse_interval(argc,argv,i));
	}
	else if (equal(argv[i],"-drift")) {
		filter.drift.push_back(parse_interval(argc,argv,i));
	}
	else if (equal(argv[i],"-nw")) {
		filter.nw = true;
	}
	else if (equal(argv[i],"-sw")) {
		filter.sw = true;
	}
    else if (equal(argv[i],"-sign_correct")) {
        filter.sign_correct = true;
    }
    else if (equal(argv[i],"-sign_correct_missing")) {
        filter.sign_correct_missing = true;
    }
	else if (equal(argv[i],"-validate")) {
		options.validate = true;
	}
	else if (equal(argv[i],"-reduction")) {
		parse_reduction_argument(argc,argv,i,options);
	}
	else if (equal(argv[i],"-stats")) {
		options.stats = true;
	}
	else if (equal(argv[i],"-model")) {
		options.model = true;
	}
	else if (equal(argv[i],"-header")) {
		action = Action::Header;
	}
	else if (equal(argv[i],"-raw")) {
		action = Action::Raw;
	}
	else if (equal(argv[i],"-raw_reduced")) {
		action = Action::RawReduced;
	}
	else if (equal(argv[i],"-reduce")) {
		action = Action::Reduce;
	}
	else if (equal(argv[i],"-test")) {
		action = Action::Test;
	}	
	else if (equal(argv[i],"-theory")) {
		parse_theory_argument(argc,argv,i,options);
	}
	else if (equal(argv[i],"-aggregate")) {
		aggregate = true;
		parse_aggregate_argument(argc,argv,i,options);
	}
	else if (equal(argv[i],"-subtract_theory")) {
		options.subtract = true;
	}
	else if (equal(argv[i],"-add_theory")) {
		options.add = true;
	}
	else if (equal(argv[i],"-data")) {
		parse_data_filename(argc,argv,i,options);
	}
	else if (equal(argv[i],"-subtract_data")) {
		options.subtract_data = true;
	}
	else if (equal(argv[i],"-subtract_model")) {
		options.subtract_model = true;
	}
	else if (equal(argv[i],"-single")) {
		options.single = true;
	}	
	else if (equal(argv[i],"-csv")) {
		options.output_format = Options::OutputFormat::CSV;
	}	
	else if (equal(argv[i],"-no_data")) {
		options.output_data = false;
	}
	else if (equal(argv[i],"-no_theory")) {
		options.output_theory = false;
	}
	else if (equal(argv[i],"-invert_data")) {
		options.invert_data = true;
	}
	else if (equal(argv[i],"-invert_theory")) {
		options.invert_theory = true;
	}
	else if (equal(argv[i],"-invert_model")) {
		options.invert_model = true;
	}
	else if (equal(argv[i],"-invert_temp")) {
		options.invert_temp = true;
	}
	else if (equal(argv[i],"-disable_desk")) {
		options.enable_desk = false;
	}
	else if (equal(argv[i],"-disable_temp")) {
		options.enable_temp = false;
	}
	else if (equal(argv[i],"-disable_earth")) {
		options.enable_earth = false;
	}
	else if (equal(argv[i],"-epoch_params")) {
		options.epoch_params = true;
	}
	else if (equal(argv[i],"-signals_dTD")) {
		parse_signals_dTD_argument(argc,argv,i,options);
	}
	else if (equal(argv[i],"-signals_ddT")) {
		parse_signals_ddT_argument(argc,argv,i,options);
	}
	else if (equal(argv[i],"-signals_dt")) {
		parse_signals_dt_argument(argc,argv,i,options);
	}
	else if (equal(argv[i],"-day_and_night")) {
		options.day_and_night = true;
	}
	else if (equal(argv[i],"-low_sun")) {
		options.low_sun = true;
	}
	else if (equal(argv[i],"-fit_amplitude")) {
		options.fit_amplitude = true;
	}
	else if (equal(argv[i],"-fit_sine")) {
		options.fit_sine = true;
	}
	else if (equal(argv[i],"-fit_disable")) {
		parse_fit_disable(argc,argv,i,options);
	}	
	else if (equal(argv[i],"-minimizer")) {
		parse_minimizer_argument(argc,argv,i,options);
	}
	else if (equal(argv[i],"-ignore")) {
		parse_ignore(argc,argv,i,options);
	}	
	else if (equal(argv[i],"-delta_chi2")) {
		parse_delta_chi_squared_argument(argc,argv,i,options);
	}
	else if (equal(argv[i],"-chi2_scale")) {
		parse_chi_squared_scale_argument(argc,argv,i,options);
	}
	else if (equal(argv[i],"-theory_params")) {
		parse_theory_params_argument(argc,argv,i,options);
	}
	else if (equal(argv[i],"-start_params")) {
		parse_start_params_argument(argc,argv,i,options);
	}
	else if (equal(argv[i],"-contour")) {
		options.contour = true;
	}
	else if (equal(argv[i],"-n")) {
		parse_index_of_refraction_argument(argc,argv,i,options);
	}
	else if (equal(argv[i],"-simulation")) {
		options.simulation = true;
	}
	else if (equal(argv[i],"-sim_seed")) {
		parse_sim_seed_argument(argc,argv,i,options);
	}
	else if (equal(argv[i],"-sim_simple")) {
		options.sim_simple = true;
	}
	else if (equal(argv[i],"-sim_sys")) {
		options.sim_sys = true;
	}
	else {
		std::cerr << "Unknown option " << argv[i] << "\n";
		throw ExitException();
	}
}



void parse_arguments(int argc,char* argv[],Filter& filter,Action& action,bool& aggregate,Options& options,
					 std::vector<std::string>& filenames)
{
	for (int i=1; i<argc; i++) {
		if (argv[i] == nullptr)
			continue;

		if (*argv[i] == '-') {// test first character
			parse_option(argc,argv,i,filter,action,aggregate,options);
		}
		else {
			filenames.push_back(argv[i]);
		}
	}
}




//-------------------------------------
// input parsing for: -aggregate fit
//-------------------------------------


IntegerInterval parse_fit_expr_epoch(const std::string& str)
{
	if (starts_with(str,"feb"))
		return {2,2};
	if (starts_with(str,"apr"))
		return {3,4};
	if (starts_with(str,"aug"))
		return {7,8};
	if (starts_with(str,"sep"))
		return {9,9};

	throw std::invalid_argument("invalid epoch");
}

IntegerInterval parse_fit_expr_interval(const std::string& str)
{
	if (!starts_with(str,"[") || !ends_with(str,"]"))
		throw std::invalid_argument("interval expected");

	auto tokens = split(str.substr(1,str.size()-2),',');
	if (tokens.size()!=2)
		throw std::invalid_argument("expected two values in an interval");

	auto from = parse_int(trim(tokens.at(0)));
	auto to = parse_int(trim(tokens.at(1)));
	if (from>to)
		throw std::invalid_argument("empty interval");

	return {from,to};
}


std::vector<IntegerInterval> parse_fit_expr_sheets(const std::string& expr)
{
	std::vector<IntegerInterval> intervals;
	auto tokens = split(expr,',','[',']');
	for (auto& token : tokens) {
		auto trimmed_token = trim(token);
		if (starts_with(trimmed_token,"[")) {
			intervals.push_back(parse_fit_expr_interval(trimmed_token));
		}
		else {
			int no = parse_int(trimmed_token);
			intervals.push_back({no,no});
		}
	}
	return intervals;
}

SignalExtractionExpression parse_fit_expression(const std::string& line)
{
	auto lr = split(line,'-');
	if (lr.size()!=2)
		throw std::invalid_argument("expected exactly one '-'");

	auto& left = lr.at(0);
	auto& right = lr.at(1);

	SignalExtractionExpression see;
	see.epoch = parse_fit_expr_epoch(left);
	see.left = parse_fit_expr_sheets(left.substr(3));
	see.right = parse_fit_expr_sheets(right);

	return see;
}



}//
