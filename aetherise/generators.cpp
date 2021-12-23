#include "generators.h"
#include "utils.h"
#include "mathematics.h"
#include "Theory.h"
#include "data_reduction.h"

#include <ctime>
#include <clocale>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>


namespace aether {







char time_of_day_symbol(TimeOfDay time_of_day)
{
	char symbol;
	switch(time_of_day) {
	case TimeOfDay::Night:
		symbol = 'n';
		break;
	case TimeOfDay::Sunrise:
	case TimeOfDay::Sunset:
		symbol = 's';
		break;
	case TimeOfDay::Day:
		symbol = 'd';
		break;
	default:
		throw std::runtime_error("unknown time of day");
	}

	return symbol;
}



std::string plusminus_symbol(const DataSheetStats& stats)
{
	std::string symbol;
	if (stats.number_of_plus > 0)
		symbol += '+';
	if (stats.number_of_minus > 0) {
		if (stats.number_of_minus > stats.number_of_plus)
			symbol = '-' + symbol;
		else
			symbol = symbol + '-';
	}
	return symbol;
}



struct __SiderealTime
{
	DataSheet::Time theta;
};



std::ostream& operator << (std::ostream& os,const __SiderealTime& st)
{
	std::stringstream ss;

	if (st.theta.hour<10)
		ss << ' ';
	ss << st.theta.hour << "ʰ";
	if (st.theta.minute<10)
		ss << ' ';
	ss << st.theta.minute << "ᵐ";

	return os << ss.str();
}



struct __SiderealTime put_sidereal_time(const DataSheet::Time& time)
{
	return {time};
}



void write_data_sheet_stats(std::ostream& os,const DataSheet& data_sheet)
{
	DataSheetStats stats = data_sheet_stats(data_sheet);
	os << data_sheet.date << ", " << std::setw(3) << data_sheet.no;
	os << std::setprecision(2) << std::fixed;
	os << "   mean T   max TD   mean dT   adjust   drift   abs drift\n";
	os << std::setw(6) << data_sheet.mean_observation_time;
	os << " " << put_sidereal_time(data_sheet.sidereal_time); // setw() does not work with varbyte unicode chars
	os << " " << time_of_day_symbol(time_of_day(data_sheet));
	os << std::setw(8) << stats.mean_T;
	os << std::setw(9) << stats.max_TD;
	os << std::setw(10) << stats.mean_dT;	
	os << std::setw(6) << plusminus_symbol(stats);
	os << std::setw(3) << stats.number_of_adjust;
	os << std::setw(8) << stats.drift;
	os << std::setw(12) << stats.abs_drift;
	os << "\n";
}



#define OUTPUT_TEMP(NAME,DIR)\
void NAME(std::ostream& os,\
		  const optional<DataSheet::Thermometers>& start,\
		  const optional<DataSheet::Thermometers>& end)\
{\
	os << std::setw(5);\
	if (start)\
		os << start->DIR;\
	else\
		os << "";\
	os << " ";\
	os << std::setw(5);\
	if (end)\
		os << end->DIR;\
	else\
		os << "";\
}

OUTPUT_TEMP(output_temp_N,N)
OUTPUT_TEMP(output_temp_E,E)
OUTPUT_TEMP(output_temp_S,S)
OUTPUT_TEMP(output_temp_W,W)



const char* location_text(const DataSheet& data_sheet)
{
	if (is_from_Cleveland(data_sheet))
		return "Cleveland";
	else 
		return "Mt. Wilson";
}



void write_header(std::ostream& os,const DataSheet& data_sheet)
{
	os << std::setprecision(1) << std::fixed;

	os << "(" << data_sheet.no << ")   " << location_text(data_sheet) << ", " << data_sheet.date << "\n";

	os << std::setw(6) << data_sheet.thermometers_start_time << " " <<
		  std::setw(5) << data_sheet.thermometers_end_time   << "\n";

	output_temp_N(os,data_sheet.thermometers_start,data_sheet.thermometers_end);

	os << "     " << std::setw(5) << data_sheet.mean_observation_time
	   << " " << std::setw(5) << data_sheet.local_mean_time << "\n";

	output_temp_E(os,data_sheet.thermometers_start,data_sheet.thermometers_end);

	os << "        θ = " << data_sheet.sidereal_time.hour << "h " << data_sheet.sidereal_time.minute << "m\n";

	output_temp_S(os,data_sheet.thermometers_start,data_sheet.thermometers_end);

	if (data_sheet.weight.has_value()) {
		os << "    weight: " << data_sheet.weight;
	}
	if (data_sheet.fringes.has_value()) {
		os << "    fringes: " << data_sheet.fringes;
	}
	os << "\n";

	output_temp_W(os,data_sheet.thermometers_start,data_sheet.thermometers_end);
	os << "\n";

	os << data_sheet.comment << "\n";

	os << data_sheet.start_time << " - " << data_sheet.end_time;
	if (data_sheet.sign_correct) {
		os << "    sign correct";
	}
	os << "\n";

	os << "\n";
}


std::string turn_attributes(const DataSheet::Turn& turn)
{
	std::string attributes;
	if (turn.adjust)
		attributes += "a";
	if (turn.plus)
		attributes += "+";
	if (turn.minus)
		attributes += "-";
	if (turn.reverse)
		attributes += "r";
	if (turn.cancel)
		attributes += "c";
	if (turn.bad)
		attributes += "b";
	if (turn.invert)
		attributes += "i";
	return attributes;
}



void write_raw(std::ostream& os,const DataSheet& data_sheet)
{
	int k=1;
	for (const DataSheet::Turn& turn : data_sheet.turns) {
		os << "# azimuth, distance\n";
		os << "0\t\"" << turn_attributes(turn) << "  " << k << "\"\n"; // columnheader
		for (int i=0;i<17;i++) {
			os << i+1 << "\t" << turn.distances.at(i) << "\n";
		}
		k++;
		write_gnuplot_data_set_separator(os);
	}
}


void write_reduce_raw(std::ostream& os,const std::vector<ReducedTurn>& reduced_turns)
{
	for (const ReducedTurn& reduced_turn : reduced_turns) {
		os << "# azimuth, distance\n";
		os << 0 << "\t" << "\"" << turn_attributes(reduced_turn.turn) << "  " << reduced_turn.i << "\"\n"; // columnheader
		for (int i=0;i<17;i++) {
			os << i+1 << "\t" << reduced_turn.displacements.at(i) << "\n";
		}

		write_gnuplot_data_set_separator(os);
	}
}



void write_reduce_stats(std::ostream& os, const DisplacementData& displacements,const Options& options)
{
	os << std::setprecision(3) << std::fixed;

	os << std::setw(15) << "";	
	os << "   mean uncertainty   χ² theory";
	if (options.model)
		os << "   χ² model ";
	os << "\n";

	os << std::setw(15) << "";
	os << std::setw(18) << mean_value(displacements.uncertainties);
	os << std::setprecision(2);
	int N = azimuths(options);
	double chi2 = chi_squared_test(displacements.data,displacements.uncertainties,displacements.theory,N);
	os << std::setw(12) << chi2;
	if (options.model) {
		chi2 = chi_squared_test(displacements.data,displacements.uncertainties,displacements.model,N);
		os << std::setw(11) << chi2;
	}
	os << "\n";
}



void write_reduced_data(std::ostream& os, const DisplacementData& displacements,
						const DataSheet& data_sheet, const Options& options)
{
	if (options.stats) {
		write_data_sheet_stats(os,data_sheet);
		write_reduce_stats(os,displacements,options);
		os << "\n";
	}
	else {
		if (options.output_format == Options::OutputFormat::CSV)
		{
			if (options.output_data) {
				os << data_sheet.date << ";" << data_sheet.no << ";";
				output_separated(os,displacements.data,";");
				os << "\n";

				os << data_sheet.date << ";" << data_sheet.no << ";";
				output_separated(os,displacements.uncertainties,";");
				os << "\n";
			}
			
			if (options.output_theory) {
				os << data_sheet.date << ";" << data_sheet.no << ";";
				output_separated(os,displacements.theory,";");
				os << "\n";
			}
			
			if (options.model) {
				os << data_sheet.date << ";" << data_sheet.no << ";";
				output_separated(os,displacements.model,";");
				os << "\n";
			}
		}
		else {
			os << "# " << data_sheet.date << ", " << data_sheet.no << "\n";
			os << "# azimuth\tdisplacement\tyerror\ttheory\tmodel\tsidereal\n";
			for (int i=0;i<17;i++) {
				os << i+1 << "\t";
				if (options.output_data) {
					os << displacements.data.at(i) << "\t";
					os << displacements.uncertainties.at(i) << "\t";
				}
				else
					os << "?\t?\t";

				if (options.output_theory)
					os << displacements.theory.at(i) << "\t";
				else
					os << "?\t";

				if (options.model)
					os << displacements.model.at(i);
				else
					os << (i==0 ? "0" : "?"); // first data zero to silence gnuplot warning
				os << "\t";

				os << time_to_h(data_sheet.sidereal_time) << "\n";
			}
			write_gnuplot_data_set_separator(os);
		}
	}
}










void write_aggregated_data(std::ostream& os, const DisplacementData& aggregated_displacements,const Options& options)
{

	if (options.stats) {		
		write_reduce_stats(os, aggregated_displacements,options);
		os << "\n";
	}
	else {
		if (options.output_format == Options::OutputFormat::CSV) {
			if (options.output_data) {
				output_separated(os,aggregated_displacements.data,";");
				os << "\n";
				output_separated(os,aggregated_displacements.uncertainties,";");
				os << "\n";
			}

			if (options.output_theory) {
				output_separated(os,aggregated_displacements.theory,";");
				os << "\n";
			}

			if (options.model) {
				output_separated(os,aggregated_displacements.model,";");
				os << "\n";
			}
		}
		else {
			os << "# azimuth\tdata\tyerror\ttheory\tmodel\n";
			for (size_t i=0; i<aggregated_displacements.data.size();i++) {
				os << i+1 << "\t";
				if (options.output_data) {
					os << aggregated_displacements.data.at(i) << "\t";
					os << aggregated_displacements.uncertainties.at(i) << "\t";
				}
				else
					os << "?\t?\t";

				if (options.output_theory)
					os << aggregated_displacements.theory.at(i) << "\t";
				else
					os << "?\t";

				if (options.model)
					os << aggregated_displacements.model.at(i);
				else
					os << (i==0 ? "0" : "?"); // first data zero to silence gnuplot warning
				os << "\n";

			}
			write_gnuplot_data_set_separator(os);
		}
	}

}







void write_aggregated_data(std::ostream& os,const std::map<int,SiderealData>& aggregated_sidereal_data,const Options& options)
{
	if (options.output_format == Options::OutputFormat::CSV) {
		using T = std::map<int,SiderealData>::const_iterator::value_type;

		output_separated(os,aggregated_sidereal_data,";",[](const T& key_value){
			return key_to_mean_sidereal_time(key_value.first);
		});
		os << "\n";

		output_separated(os,aggregated_sidereal_data,";",[](const T& key_value){
			return key_value.second.data_amplitude;
		});
		os << "\n";

		output_separated(os,aggregated_sidereal_data,";",[](const T& key_value){
			return key_value.second.theory_amplitude;
		});
		os << "\n";

		if (options.model) {
			output_separated(os,aggregated_sidereal_data,";",[](const T& key_value){
				return key_value.second.model_amplitude;
			});
			os << "\n";
		}
	}
	else {
		os << "# sidereal\tdata\ttheory\tmodel\tTD\n";

		double last_sidereal = 0;
		if (!aggregated_sidereal_data.empty()) {
			last_sidereal = key_to_mean_sidereal_time(aggregated_sidereal_data.begin()->first);
		}
		int i=0;
		for (const auto& key_value : aggregated_sidereal_data) {
			auto sidereal_time = key_to_mean_sidereal_time(key_value.first);
			if (sidereal_time-last_sidereal > 1)
				os << "\n"; // don't connect these points with lines

			os << sidereal_time << "\t";
			os << key_value.second.data_amplitude << "\t";
			os << key_value.second.theory_amplitude << "\t";
			if (options.model)
				os << key_value.second.model_amplitude << "\t";
			else
				os << (i==0 ? "0" : "?") << "\t"; // first data zero to silence gnuplot warning
			os << key_value.second.TD*10 << "\n";

			last_sidereal = sidereal_time;
			i++;
		}
		write_gnuplot_data_set_separator(os); // usefull to concat files
	}
}




void write_list(std::ostream& os,const std::vector<DataSheet>& data_sheets,const Options& options)
{
	os << "\"date\";\"no\";\"time\";\"sidereal time\";\"time of day\";\"duration\";\"turns\";\"duration/turn\";\"attributes\";\"weather\";\"uncertainty\";\"T\";\"TD\";\"dT\";"
		  "\"sign correct\";\"adjust\";\"sign\";\"drift\";\"abs drift\";\"notes\"\n";

	for (const DataSheet& data_sheet : data_sheets) {
		auto reduced_data = reduce_data(data_sheet,options);
		auto stats = data_sheet_stats(data_sheet);
		auto mean_u = mean_value(reduced_data.uncertainties);
		auto duration = periodic_distance(time_to_h(data_sheet.start_time),time_to_h(data_sheet.end_time),24) * 60; // min
		auto duration_per_turn = duration/(data_sheet.turns.size()+stats.number_of_adjust); // min
		
		
		output_args_separated(os,";",
							  data_sheet.date,data_sheet.no,
							  data_sheet.mean_observation_time,data_sheet.sidereal_time,
							  time_of_day_symbol(time_of_day(data_sheet)),
							  duration,
							  data_sheet.turns.size(),
							  duration_per_turn,							  
							  quote(data_sheet_attributes(data_sheet)),
							  quote(data_sheet_weather(data_sheet)),
							  mean_u,stats.mean_T,stats.max_TD,stats.mean_dT,
							  (data_sheet.sign_correct ? "x" : ""),
							  stats.number_of_adjust,quote(plusminus_symbol(stats)),
							  stats.drift,stats.abs_drift,
							  quote(data_sheet.comment));
		os << "\n";
	}
}



void write_chi_squared_stats(std::ostream& os,double chi,int f,double p,const std::string& prefix)
{
	os << std::setprecision(3) << std::fixed;
	os << prefix << "χ² = " << chi << "\n";
	os << prefix << "f  = " << f << "\n";
	os << prefix << "χ²/f = " << chi/f << "\n";

	os << std::setprecision(6) << std::defaultfloat;
	os << prefix << "p-value = " << p << "\n";
}



void write_spectrum(std::ostream& os,const std::map<double,std::complex<double>>& results)
{	
	os << "# freq\tamplitude\n";
	for (auto& entry : results) { // map is sorted		
		auto& f = entry.first;
		auto& z = entry.second;
		 os << f << "\t" << std::abs(z)/20 << "\n";		 
	}
	
	write_gnuplot_data_set_separator(os);
}

}//
