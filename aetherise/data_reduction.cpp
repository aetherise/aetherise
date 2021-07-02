
#include "data_reduction.h"
#include "utils.h"
#include "mathematics.h"
#include "mini.h"

#include <algorithm>
#include <numeric>
#include <iostream>
#include <random>


namespace aether {



std::array<double,17> reduce_from_drift_and_offset(const std::array<float,17>& distances)
{
	std::array<double,17> displacements;

	// remove (assumed) linear drift
	auto delta = distances.back()-distances.front();
	for(size_t i=0;i<17;i++) {
		displacements[i] = distances[i] - delta/16.0*i;
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
 * @param single reduce separate turns to single period before calculating uncertainties?
 * @return
 */
std::array<double,17>
standard_uncertainties(const DataSheet& data_sheet,const Options& options)
{
	std::vector<std::array<double,17>> turns;
	std::array<double,17> estimates {};

	int n = 0;
	selected_and_transformed_turns(data_sheet,options,[&](int,const DataSheet::Turn&,const std::array<float,17>& distances) {
		std::array<double,17> displacements = reduce_from_drift_and_offset(distances);

		if (options.single)
			reduce_to_single_period(displacements);

		add_array(estimates,displacements);
		turns.push_back(std::move(displacements));
		n++;
	});

	for (auto& x : estimates)
		x /= n;


	// sample variance
	std::array<double,17> s {};
	for (auto& displacements : turns) {
		for (size_t i=0;i<17;i++) {
			s.at(i) += sqr(displacements.at(i)-estimates.at(i));
		}
	}	
	for (auto& si : s)
		si /= (n-1);

	// standard uncertainty
	for (auto& si : s)
		si = std::sqrt(si/n);

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







ReducedData MillersReduction::reduce(const DataSheet &data_sheet, const Options &options) const
{
	ReducedData reduced_data {};

	int used_turns = 0;
	// step 1: column sum
	selected_and_transformed_turns(data_sheet,options,[&](int,const DataSheet::Turn&,const std::array<float,17>& distances){
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
	
	if (options.single)
		reduce_to_single_period(reduced_data.displacements); // Miller did this after step 4, but its better before
	
	reduced_data.uncertainties = standard_uncertainties(data_sheet,options);

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
		
	selected_and_transformed_turns(data_sheet,options,[&](int,const DataSheet::Turn&,const std::array<float,17>& distances){
		
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


	reduced_data.uncertainties = standard_uncertainties(data_sheet,options);

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



Estimate<std::complex<double>> reduce_ks(const std::vector<std::complex<double>>& ks)
{
	auto R = estimate_expected_value(ks.begin(),ks.end(),[](const std::complex<double>& z){
		return z.real();
	});	
	
	auto I = estimate_expected_value(ks.begin(),ks.end(),[](const std::complex<double>& z){
		return z.imag();
	});	
		
	std::complex<double> z {R.m,I.m};
	std::complex<double> u {R.u,I.u};
	
	return {z,u};		
}



ReducedData DFTReduction::reduce(const DataSheet &data_sheet, const Options &options) const
{
	ReducedData reduced_data {};
				
	std::vector<std::complex<double>> k1s;
	std::vector<std::complex<double>> k2s;
			
	selected_and_transformed_turns(data_sheet,options,[&](int,const DataSheet::Turn&,const std::array<float,17>& distances){		
		// step 1: remove (assumed) linear drift
		std::array<double,17> displacements = reduce_from_drift_and_offset(distances);
		for (auto& d : displacements)
			d /= 20.; // to wave length
		if (options.single)
			reduce_to_single_period(displacements);

		auto k1 = DFT_analyze(1,displacements.begin(),displacements.end()-1);
		auto k2 = DFT_analyze(2,displacements.begin(),displacements.end()-1);		
						
		k1s.push_back(k1);
		k2s.push_back(k2);				
	});
		
	auto z1 = reduce_ks(k1s);
	auto z2 = reduce_ks(k2s);
				
	reduced_data.z1 = z1;
	reduced_data.z2 = z2;	
	
			
	// generate data for harmonics 1 and 2			
	const double h = 1e-6;
		
	for (size_t i=0;i<17;i++) {		
		auto y1 = propagate_uncertainties(Real(z1),Imag(z1),[&i](double R,double I){
			auto z = std::complex<double>{R,I};
			auto A = std::abs(z);
			auto phi = std::arg(z);
			return A*std::cos(i*AETHER_2PI/16 + phi);
		},h);
				
		auto y2 = propagate_uncertainties(Real(z2),Imag(z2),[&i](double R,double I){
			auto z = std::complex<double>{R,I};
			auto A = std::abs(z);
			auto phi = std::arg(z);
			return A*std::cos(i*AETHER_2PI/8 + phi);
		},h);
		
		reduced_data.displacements[i] = y1.m + y2.m;
		reduced_data.uncertainties[i] = std::sqrt(sqr(y1.u)+sqr(y2.u));
	}
	
	return reduced_data;
}



ReducedData Roberts2006::reduce(const DataSheet &data_sheet, const Options &options) const
{
	ReducedData reduced_data {};

	
	std::vector<std::array<float,17>> turns;
	std::vector<float> model; // helper, "turns" as a sequence
			
	selected_and_transformed_turns(data_sheet,options,[&](int,const DataSheet::Turn&,const std::array<float,17>& distances){
		turns.push_back(distances);
	});

	const auto first_turn = turns.at(0); 
		
	// subtract first half turn 
	if (options.subtract_first_half_turn) {
		std::cerr << "First half turn will be subtracted\n";		
		for (auto& turn : turns) {
			for (size_t i=0;i<8;i++) {
				turn.at(i)   -= first_turn.at(i);
				turn.at(i+8) -= first_turn.at(i);
			}
			turn.at(16) -= first_turn.at(0);
		}
	}
		
	// de-adjust
	float offset = 0;
	for (auto& turn : turns) {
		auto t0 = turn.at(0);
		for (auto& ti : turn)
			ti -= t0 - offset; // de-adjust, remove jumps		
		offset += turn[16]-turn[0];
							
		for (auto it = turn.begin();it!=turn.end()-1;++it) 
			model.push_back(*it);			
	};
	
	// check stability
	size_t unstable_count = 0;
	for (const auto& turn : turns) {
		for (size_t i=0;i<9;i++) {
			if (std::abs(turn.at(i+8)-turn.at(i)) >= 5) {// half fringe limit
				unstable_count++;
				break;
			}
		}
	}
	std::cerr << "Number of stable turns: " << turns.size()-unstable_count << " out of " << turns.size() << "\n";
	
	
	// Chi squared function
	auto chi2f = [&model](const std::array<float,7>& x){
		float chi2 = 0;
		for (size_t i=1;i<model.size();i++) {
			auto p0 = (i-1)%8==0 ? 0 : x.at((i-1)%8-1);
			auto p1 = (i)%8==0 ? 0 : x.at((i)%8-1);
			chi2 += sqr((model.at(i)+p1)-(model.at(i-1)+p0));
		}
		return chi2;
	};
	

	// find minimum simply by permutating the 7 parameters
	std::array<float,7> params {};
	auto min_params = params;
	float min_chi2 = INFINITY;
	
	float p_delta = 4;
	float pfrom = -p_delta;
	float pto = p_delta;
	for (int i=0;i<4;i++) {
		float prev_min_chi = min_chi2;
		std::cerr << "Iteration " << i+1 << std::flush;
		for (auto p1=pfrom;p1<=pto;p1++) {
			std::cerr << ".";
			for (auto p2=pfrom;p2<=pto;p2++) // yes, it's a for-loop with float			
			for (auto p3=pfrom;p3<=pto;p3++)
			for (auto p4=pfrom;p4<=pto;p4++)
			for (auto p5=pfrom;p5<=pto;p5++)
			for (auto p6=pfrom;p6<=pto;p6++)
				for (auto p7=pfrom;p7<=pto;p7++) {
					std::array<float,7> ps {p1,p2,p3,p4,p5,p6,p7};
					add_array(ps,params);
					auto chi2 = chi2f(ps);
					if (chi2<min_chi2) {
						min_chi2 = chi2;
						min_params = ps;
					}	
				}
		}
		std::cerr << "\n";
		
		params = min_params; // new starting point for next iteration	
		if (std::abs(prev_min_chi - min_chi2) < 0.1f)
			break;
	}
	
	
	
	std::cerr << "Model parameters: ";
	output_separated(std::cerr,params,", ");
	std::cerr << "\n";
	std::cerr << "χ² = " << min_chi2 << "\n";
	
	
	// Create the model with the found parameters and do the DFT
	DFTGoertzel dft({2},16);
	for (auto& turn : turns) {
		for (size_t i=0;i<7;i++) {
			turn.at(1+i) += params.at(i);
			turn.at(9+i) += params.at(i);
		}	
				
		dft.analyze(turn.begin(),turn.end()-1);			
	}
	
	auto result = dft.result();
	std::cerr << "1/2 turn DFT amplitude = " << std::abs(result.at(2))/20 << "\n"; // 1/10 fringe to wave length
	
	
	
	//-----------------------------------------------------------
	// reduction - with Roberts' model
	//-----------------------------------------------------------
	
	int used_turns = 0;
	for (auto& turn : turns) {					
		// step 1: remove (assumed) linear drift
		std::array<double,17> displacements = reduce_from_drift_and_offset(turn);
		if (options.single)
			reduce_to_single_period(displacements);
				
		// step 2: add turn
		std::transform(reduced_data.displacements.begin(),reduced_data.displacements.end(),displacements.begin(),
					   reduced_data.displacements.begin(), std::plus<double>());
		used_turns++;
	};
	

	// step 3: divide by number of turns/rows (usually 20)
	for (double& d : reduced_data.displacements) {
		d /= used_turns;
	}


	reduced_data.uncertainties = standard_uncertainties(data_sheet,options);

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




std::unique_ptr<ReductionMethod> create_reduction_method(Options::ReductionMethod method)
{
	std::unique_ptr<ReductionMethod> reduction_method;

	switch (method) {
	case Options::ReductionMethod::Miller:
		reduction_method = make_unique<MillersReduction>();
		break;		
	case Options::ReductionMethod::DFT:
		reduction_method = make_unique<DFTReduction>();
		break;
	case Options::ReductionMethod::Roberts2006:
		reduction_method = make_unique<Roberts2006>();
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
