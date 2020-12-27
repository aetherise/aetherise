
#include "data_reduction.h"
#include "utils.h"
#include "mathematics.h"

#include <algorithm>
#include <numeric>
#include <iostream>

namespace aether {



std::array<double,17> reduce_from_drift_and_offset(const std::array<short,17>& distances)
{
	std::array<double,17> displacements;

	// remove (assumed) linear drift
	double delta = distances.back()-distances.front();
	for(size_t i=0;i<17;i++) {
		displacements.at(i) = distances.at(i) - delta/16.0*i;
	}

	// move by mean ordinate
	double mean_ordinate = std::accumulate(displacements.begin(),displacements.end()-1,0.0) / 16.0;
	for (double& d : displacements) {
		d -= mean_ordinate;
	}

	return displacements;
}




/**
 * \~german
 * Standardmessunsicherheiten
 *
 * \~english
 *
 * \~
 * @param data_sheet
 * @param options
 * @param single reduce separate turns to single period before calculating uncertainies?
 * @return
 */
std::array<double,17>
standard_uncertainties(const DataSheet& data_sheet,const Options& options,bool single)
{
	std::vector<std::array<double,17>> turns;
	std::array<double,17> estimates {};

	int n = 0;
	selected_and_transformed_turns(data_sheet,options,[&](int,const DataSheet::Turn&,const std::array<short int,17>& distances) {
		std::array<double,17> displacements = reduce_from_drift_and_offset(distances);

		if (single)
			reduce_to_single_period(displacements);

		add_array(estimates,displacements);
		turns.push_back(std::move(displacements));
		n++;
	});

	for (auto& x : estimates)
		x /= n;


	std::array<double,17> s {};
	for (auto& displacements : turns) {
		for (size_t i=0;i<17;i++) {
			s.at(i) += sqr(displacements.at(i)-estimates.at(i));
		}
	}

	// sample (experimental) standard deviation
	for (auto& si : s)
		si = std::sqrt(si/(n-1));

	// standard uncertainty
	for (auto& si : s)
		si /= std::sqrt(n);

	if (single) {
		// Correction of uncertainty, because we evaluate the uncertainty
		// with values which uncertainties are lowered (approximately) by the factor sqrt(2)/2.
		// TODO Unsure if this is correct.
		for (auto& si : s)
			si *= 2./std::sqrt(2);
	}

	for (auto& si : s)
		si *= options.chi_squared_scale;
	
	return s;
}



void reduce_to_single_period(std::array<double,17>& displacements)
{
	std::transform(displacements.begin(),displacements.begin()+8,displacements.begin()+8,
				   displacements.begin(),[](double d1,double d2){
		return (d1+d2)*0.5;
	});
	displacements.at(8) = (displacements.at(8) + displacements.at(16))*0.5;
	std::copy(displacements.begin()+1,displacements.begin()+9,displacements.begin()+9);
}



void reduce_to_single_period(std::array<double,17>& displacements,std::array<double,17>& uncertainties)
{
	reduce_to_single_period(displacements);

	std::transform(uncertainties.begin(),uncertainties.begin()+8,uncertainties.begin()+8,
				   uncertainties.begin(),[](double u1,double u2){
		return std::sqrt(sqr(u1)+sqr(u2))*0.5;
	});
	uncertainties.at(8) = std::sqrt(sqr(uncertainties.at(8)) + sqr(uncertainties.at(16)))*0.5;
	std::copy(uncertainties.begin()+1,uncertainties.begin()+9,uncertainties.begin()+9);
}



void reduce_to_single_period(ReducedData& reduced_data)
{
	reduce_to_single_period(reduced_data.displacements, reduced_data.uncertainties);
}



ReducedData MillersReduction::reduce(const DataSheet &data_sheet, const Options &options) const
{
	ReducedData reduced_data {};

	int used_turns = 0;
	// step 1: column sum
	selected_and_transformed_turns(data_sheet,options,[&](int,const DataSheet::Turn&,const std::array<short int,17>& distances){
		std::transform(reduced_data.displacements.begin(),reduced_data.displacements.end(),distances.begin(),
					   reduced_data.displacements.begin(),std::plus<double>());
		used_turns++;
	});

	// step 2: remove (assumed) linear drift
	double delta = reduced_data.displacements.back()-reduced_data.displacements.front();
	for(size_t i=1;i<17;i++) {
		reduced_data.displacements.at(i) -= delta/16.0*i;
	}
	// step 3: divide by number of turns/rows (usually 20)
	for (double& d : reduced_data.displacements) {
		d /= used_turns;
	}

	reduced_data.uncertainties = standard_uncertainties(data_sheet,options,false);
	if (options.single)
		reduce_to_single_period(reduced_data); // Miller did this after step 4, but its better before

	// step 4: move by mean ordinate
	double mean_ordinate = std::accumulate(reduced_data.displacements.begin(),reduced_data.displacements.end()-1,0.0) / 16.0;
	for (double& d : reduced_data.displacements) {
		d -= mean_ordinate;
	}	

	// step 5: in wave length
	for (double& d : reduced_data.displacements) {
		d /= 20.0;
	}

	for (double& u : reduced_data.uncertainties) {
		u /= 20.0;
	}


	return reduced_data;
}





ReducedData SeparateReduction::reduce(const DataSheet &data_sheet, const Options &options) const
{
	ReducedData reduced_data {};

	int used_turns = 0;

	selected_and_transformed_turns(data_sheet,options,[&](int,const DataSheet::Turn&,const std::array<short int,17>& distances){
		// step 1: remove (assumed) linear drift
		std::array<double,17> displacements = reduce_from_drift_and_offset(distances);
		if (options.single)
			reduce_to_single_period(displacements);

		// step 2: add turn
		std::transform(reduced_data.displacements.begin(),reduced_data.displacements.end(),displacements.begin(),
					   reduced_data.displacements.begin(), std::plus<double>());
		used_turns++;
	});


	// step 3: divide by number of turns/rows (usually 20)
	for (double& d : reduced_data.displacements) {
		d /= used_turns;
	}


	reduced_data.uncertainties = standard_uncertainties(data_sheet,options,options.single);

	// step 4: move by mean ordinate
	double mean_ordinate = std::accumulate(reduced_data.displacements.begin(),reduced_data.displacements.end()-1,0.0) / 16.0;
	for (double& d : reduced_data.displacements) {
		d -= mean_ordinate;
	}

	// step 5: in wave length
	for (double& d : reduced_data.displacements) {
		d /= 20.0;
	}

	for (double& u : reduced_data.uncertainties)
		u /= 20.0;

	return reduced_data;
}




std::unique_ptr<ReductionMethod> create_reduction_method(Options::DataReductionMethod method)
{
	std::unique_ptr<ReductionMethod> reduction_method;

	switch (method) {
	case Options::DataReductionMethod::Miller:
		reduction_method = make_unique<MillersReduction>();
		break;	
	case Options::DataReductionMethod::Separate:
		reduction_method = make_unique<SeparateReduction>();
		break;
	default:
		throw std::runtime_error("unknown reduction method");
	}

	return reduction_method;
}



ReducedData reduce_data(const DataSheet& data_sheet,const Options& options)
{
	auto reduction_method = create_reduction_method(options.reduction_method);

	return reduction_method->reduce(data_sheet,options);
}


}//aether
