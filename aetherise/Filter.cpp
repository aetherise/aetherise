
#include "Filter.h"
#include "data_reduction.h"
#include "Theory.h"

#include <cmath>

namespace aether {

// prototype
std::array<double,17> fringe_displacements(const Theory& theory,const TheoryParameters& params,
										   const DataSheet& data_sheet,const Options& options);


bool is_filter_set(const Filter::Interval& range)
{
	return !std::isinf(range.min) || !std::isinf(range.max);
}



bool is_filter_set(const std::vector<Filter::Interval>& ranges)
{
	for (const auto& range : ranges)
		if (is_filter_set(range))
			return true;
	return false;
}



struct not_selected {}; // used as exception

template<typename T>
void check_selected(T value, const std::vector<Filter::Interval>& ranges)
{
	if (ranges.empty())
		return;
	for (const auto& range : ranges) {
		if (value >= range.min && value <= range.max)
			return;
	}

	throw not_selected();
}



template<>
void check_selected(DataSheet::Thermometers value, const std::vector<Filter::Interval>& ranges)
{
	check_selected(value.N,ranges);
	check_selected(value.E,ranges);
	check_selected(value.S,ranges);
	check_selected(value.W,ranges);
}



bool selected(const DataSheet& data_sheet,const Options& options,const Filter& filter)
{
	try {
		//--------------
		// sheet values
		//--------------

		check_selected(data_sheet.no,filter.no);

		check_selected(data_sheet.date.year,filter.year);
		check_selected(data_sheet.date.month,filter.month);
		check_selected(data_sheet.date.day,filter.day);

		check_selected(time_to_h(data_sheet.mean_observation_time),filter.time);
		check_selected(time_to_h(data_sheet.sidereal_time),filter.sidereal_time);



		if (is_filter_set(filter.weight)) {
			if (data_sheet.weight.has_value()) {
				check_selected(*data_sheet.weight,filter.weight);
			}
			else return false;
		}

		if (is_filter_set(filter.fringes)) {
			if (data_sheet.fringes.has_value()) {
				check_selected(*data_sheet.fringes,filter.fringes);
			}
			else return false;
		}



		// T
		if (is_filter_set(filter.T)) {
			if (data_sheet.thermometers_start || data_sheet.thermometers_end) {
				if (data_sheet.thermometers_start) {
					check_selected(*data_sheet.thermometers_start,filter.T);
				}
				if (data_sheet.thermometers_end) {
					check_selected(*data_sheet.thermometers_end,filter.T);
				}
			}
			else
				return false;

		}

		// TD
		if (is_filter_set(filter.TD)) {
			if (data_sheet.thermometers_start || data_sheet.thermometers_end) {
				if (data_sheet.thermometers_start) {
					auto D = max_T(*data_sheet.thermometers_start) - min_T(*data_sheet.thermometers_start);
					check_selected(D,filter.TD);
				}
				if (data_sheet.thermometers_end) {
					auto D = max_T(*data_sheet.thermometers_end) - min_T(*data_sheet.thermometers_end);
					check_selected(D,filter.TD);
				}
			}
			else
				return false;
		}

		//dT
		if (is_filter_set(filter.dT)) {
			if (!(data_sheet.thermometers_start && data_sheet.thermometers_end)) {
				return false;
				// after this, temperature value access is allowed
			}

			check_selected(data_sheet.thermometers_end->N - data_sheet.thermometers_start->N,filter.dT);
			check_selected(data_sheet.thermometers_end->E - data_sheet.thermometers_start->E,filter.dT);
			check_selected(data_sheet.thermometers_end->S - data_sheet.thermometers_start->S,filter.dT);
			check_selected(data_sheet.thermometers_end->W - data_sheet.thermometers_start->W,filter.dT);
		}



		//-----------
		// stats
		//-----------

		DataSheetStats stats = data_sheet_stats(data_sheet);
		check_selected(stats.drift,filter.drift);
		check_selected(stats.number_of_adjust,filter.adjust);
		if (is_filter_set(filter.mean_dT)) {
			if (stats.mean_dT.has_value())	{
				check_selected(*stats.mean_dT,filter.mean_dT);
			}
			else
				return false;
		}

		//---------------
		// computed
		//---------------

		if (is_filter_set(filter.amplitude) || is_filter_set(filter.uncertainty)) {
			ReducedData reduced_data = reduce_data(data_sheet,options);

			if (is_filter_set(filter.amplitude)) {				
				double amplitude = max_abs_value(reduced_data.displacements);
				check_selected(amplitude,filter.amplitude);
			}

			if (is_filter_set(filter.uncertainty)) {
				auto mean = mean_value(reduced_data.uncertainties);
				check_selected(mean,filter.uncertainty);
			}
		}

		if (is_filter_set(filter.theory_amplitude)) {
			auto theory = create_theory(options.theory);
			auto displacements = fringe_displacements(*theory,CMB_dipole,data_sheet,options);
			double amplitude = max_abs_value(displacements);

			check_selected(amplitude,filter.theory_amplitude);
		}


		//---------
		// options
		//---------
		if (filter.nw && data_sheet.desk_in_sw)
			return false;
		if (filter.sw && !data_sheet.desk_in_sw)
			return false;
        if (filter.sign_correct && !data_sheet.sign_correct)
            return false;
        if (filter.sign_correct_missing && data_sheet.sign_correct)
            return false;


	}
	catch(not_selected&) {
		return false;
	}

	return true;
}



}//aether
