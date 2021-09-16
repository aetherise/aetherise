
#include "aetherise.h"

#include "DataSheet.h"
#include "utils.h"
#include "Filter.h"
#include "data_reduction.h"
#include "Theory.h"
#include "astro.h"
#include "models.h"
#include "generators.h"
#include "cmd_line.h"
#include "mathematics.h"
#include "mini.h"
#include "physics.h"

#include <iostream>
#include <cmath>
#include <memory>
#include <algorithm>
#include <map>
#include <iomanip>
#include <queue>
#include <functional>
#include <random>
#include <sstream>


namespace aether {



/**
 * \~german
 * Die mit -data angebenen Daten von den gegebenen Datensätzen subtrahieren.
 * 
 * \~english
 * Subtract the data given with -data from the given data records.
 * 
 * \~
 * @param displacements
 * @param options
 */
void subtract_data(DisplacementData& displacements,const Options& options)
{
	// TODO: data rows should be generic? No handling for -no_data currently!
	if (options.data1.size()<17 || options.data2.size()<17 
			|| (options.output_theory && options.data3.size()<17) 
			|| (options.model && options.data4.size()<17)) {		
		std::cerr << "ERROR: " << (options.model ? 4 : 3) << " data rows with 17 values needed\n";
		throw ExitException();
	}

	std::transform(displacements.data.begin(),displacements.data.end(),options.data1.begin(),
				   displacements.data.begin(),std::minus<double>());

	std::transform(displacements.uncertainties.begin(),displacements.uncertainties.end(),options.data2.begin(),
				   displacements.uncertainties.begin(),[](double u1,double u2){
			return std::sqrt(sqr(u1) + sqr(u2));
		});

	if (options.output_theory) {
		std::transform(displacements.theory.begin(),displacements.theory.end(),options.data3.begin(),
					   displacements.theory.begin(),std::minus<double>());
	}

	if (options.model) {
		std::transform(displacements.model.begin(),displacements.model.end(),options.data4.begin(),
					   displacements.model.begin(),std::minus<double>());
	}
}



/**
 * \~german
 * Erdbahn addieren
 *
 * Addiert den Geschwindigkeitsvektor der Erde, auf ihrer Bahn um die Sonne,
 * zu dem gegebenen Geschwindigkeitsvektor des Sonnensystems.
 *
 * \~english
 * Adds the velocity vector of the earth, on its orbit around the sun,
 * to the given velocity vector of the solar system.
 *
 * \~
 * @param params
 * @param data_sheet
 * @return modified parameters
 */
TheoryParameters add_earth_orbit(const TheoryParameters& params,const DataSheet& data_sheet,const Options& options)
{
	Calendar cal = calendar_date(data_sheet);
	auto JD = julian_date(cal);
	JD -= rad_to_h(options.longitude)/24; // to UT, add longitude as time (makes almost no difference)			
	Equatorial apex = earth_apex(JD);

	Polar3 earth;
	earth.r = Earth_average_orbital_speed;
	earth.th = rad(90)-apex.de;
	earth.ph = apex.ra;

	Polar3 sun;
	sun.r = params.v;
	sun.th = rad(90)-params.d;
	sun.ph = params.a;

	auto v = vector3(sun) + vector3(earth);
	auto p = polar3(v);

	TheoryParameters r;
	r.v = p.r;
	r.a = p.ph;
	r.d = rad(90)-p.th;
	return r;
}



/**
 * \~german
 * Luftfeuchtigkeit
 *
 * Relative Luftfeuchtigkeit aus den Wetterdaten eines Datenblattes abschätzen.
 *
 * \~english
 * Estimate the relative humidity using the weather data of the given data sheet.
 *
 * \~
 * @param data_sheet
 * @return A value in the interval [0.0, 1.0]
 */
double humidity(const DataSheet& data_sheet)
{
	double h;

	// all values are just estimates
	switch(data_sheet.weather_obscuration)
	{
	case DataSheet::WeatherObscurations::Clear:
		h = 0.5;
		break;
	case DataSheet::WeatherObscurations::Haze:
		h = 0.5; // TODO: what exactly is haze!?
		break;
	case DataSheet::WeatherObscurations::Mist:
		h = 0.9;
		break;
	case DataSheet::WeatherObscurations::Fog:
		h = 1.0;
		break;
	default:
		throw std::runtime_error("unknown atmospheric obscuration");
	}

	if (data_sheet.weather_rain) {
		if (h < 0.8)
			h = 0.8;
	}

	return h;
}



/**
 * \~german
 * Brechungsindex
 *
 * Den Brechungsindex von Luft auf dem Mount Wilson aus den Wetterdaten
 * des gegebenen Datenblattes berechnen.
 *
 * \~english
 * Calculates the refractive index of air on top of Mount Wilson
 * using the weather data of the given data sheet.
 *
 * \~
 * @param data_sheet
 * @return
 */
double index_of_refraction(const DataSheet& data_sheet)
{
	// The CO2 value (ppm) of the year 1925 is from:
	// D.M.Etheridge et al.: Historical CO2 Records from the Law Dome DE08, DE08-2, and DSS Ice Cores (1006 A.D.-1978 A.D).
	// Carbon Dioxide Information Analysis Center, 1998
	// doi:10.3334/CDIAC/ATG.011
	// https://cdiac.ess-dive.lbl.gov/ftp/trends/co2/lawdome.smoothed.yr20

	auto t = mean_T(data_sheet);
	if (!t)
		return 1.00023; // estimated mean for the best groups in September and February

	auto T = in_Kelvin(*t);
	// refractive index of air in 1700 m altitude
	auto p = barometric_formula(1700,101325,T+11.5); // 11.5=1700*0.0065 (temperature gradient)
	auto lambda = Millers_Interferometer_Wave_Length*1e+6; // (µm)
	return refractive_index_of_air(p,T,lambda,humidity(data_sheet),305);
}



/**
 * \~german Streifenverschiebung
 * Streifenverschiebung berechnen unter Beachtung verschiedener Einstellungen.
 * 
 * \~english
 * Calculate the fringe displacements considering different settings.
 * 
 * \~
 * @param theory
 * @param params
 * @param data_sheet
 * @param options
 * @return 
 */
std::array<double,17> fringe_displacements(const Theory& theory,const TheoryParameters& params,
										   const DataSheet& data_sheet,const Options& options)
{	
	auto tparams = options.enable_earth ? add_earth_orbit(params,data_sheet,options) : params;

	double n;
	if (std::isnan(options.index_of_refraction))
		n = index_of_refraction(data_sheet);
	else
		n = options.index_of_refraction;

	auto theta = h_to_rad(time_to_h(data_sheet.sidereal_time));	
	//auto theta = sidereal_time(calendar_date(data_sheet),0);
	//theta += rad(8*15)+options.longitude; // correction for timezone and longitude
	bool invert = options.invert_theory != data_sheet.inverted;	
	return fringe_displacements(theory,tparams,options.latitude,n,theta,invert);
}



/**
 * \~german Datenblatt verarbeiten
 * Führt die Datenreduzierung und weitere über die Schalter ausgewählte Datenverarbeitung aus.
 * Erzeugt auch das erwartete theoretische Signal und das Modell.
 * 
 * \~english
 * Executes the data reduction and further data processing selected by the options.
 * Also generates the theoretical signal und the model.
 *  
 * \~
 * @param data_sheet
 * @param options
 * @return 
 */
DisplacementData process_data_sheet(const DataSheet& data_sheet, const Options& options)
{
	DisplacementData displacements;

	const ReducedData reduced_data = reduce_data(data_sheet,options);
	displacements.data = reduced_data.displacements;
	displacements.uncertainties = reduced_data.uncertainties;	
	displacements.z1 = reduced_data.z1;
	displacements.z2 = reduced_data.z2;

	std::unique_ptr<Theory> theory = create_theory(options.theory);	
	displacements.theory = fringe_displacements(*theory,options.theory_params,data_sheet,options);

	if (options.model) {
		displacements.model = systematic_error_displacements(data_sheet,options);
	}
	else {
		displacements.model = {};
	}

	if (options.subtract_model) {
		sub_array(displacements.data,displacements.model);
		sub_array(displacements.model,displacements.model);
	}

	add_array(displacements.model,displacements.theory);


	if (options.add) {
		add_array(displacements.data,displacements.theory);
		//add_array(displacements.model,displacements.theory); // makes no sense
	}

	if (options.subtract) {
		sub_array(displacements.data,displacements.theory);
		sub_array(displacements.model,displacements.theory);
	}

	if (options.subtract_data)
		subtract_data(displacements,options);



	return displacements;
}




void execute_reduce(const DataSheet& data_sheet,const Options& options)
{	
	DisplacementData displacements = process_data_sheet(data_sheet,options);

	write_reduced_data(std::cout,displacements,data_sheet,options);
}




std::string n_of_m(int n,int m,bool space=true)
{
	std::ostringstream oss;
	oss << n << (space ? " / " : "/") << m;
	return oss.str();
}




size_t test_for_normality(const DataSheet& data_sheet,const Options& options, 
						  std::vector<double>& As,std::vector<int>& n_sig)
{
	const auto N = azimuths(options);
	std::vector<std::vector<double>> azimuth_samples(N);

	selected_and_transformed_turns(data_sheet,options,[&](int ,const DataSheet::Turn& ,const std::array<float,17>& distances){
		auto displacements = reduce_from_drift_and_offset(distances);

		if (options.single)
			reduce_to_single_period(displacements);

		for (int i=0;i<N;i++) {
			azimuth_samples.at(i).push_back(displacements.at(i));
		}
	});
	const auto n = azimuth_samples.at(0).size(); // same for all


	// Anderson-Darling Test for each azimuth.
	// D'Agostino is closer to Shapiro-Wilk test results (in R 3.3.3) than Stephens.
	auto type = ADTestType::DAgostino;
	for (auto& samples : azimuth_samples) {
		std::sort(samples.begin(),samples.end());
		As.push_back( test_for_normality(samples,type) );
	}

	// count rejected
	int n_50 = 0;
	int n_25 = 0;
	int n_15 = 0;
	int n_10 = 0;
	int n_5 = 0;
	int n_1 = 0;

	auto quantiles = test_quantiles(type);
	for (const auto& A : As) {
		if (A > quantiles.q_1) {
			n_50++; n_25++; n_15++; n_10++; n_5++; n_1++;
		}
		else if (A > quantiles.q_5) {
			n_50++; n_25++; n_15++; n_10++; n_5++;
		}
		else if (A > quantiles.q_10) {
			n_50++; n_25++; n_15++; n_10++;
		}		
		else if (A > quantiles.q_15) {
			n_50++; n_25++; n_15++;
		}
		else if (A > quantiles.q_25) {
			n_50++; n_25++;
		}
		else if (A > quantiles.q_50) {
			n_50++;
		}
	}

	n_sig = {n_50, n_25, n_15, n_10, n_5, n_1};
	return n;
}




void execute_test(const DataSheet& data_sheet,const Options& options)
{
	const auto N = azimuths(options);

	std::vector<double> As;
	std::vector<int> n_sig;
	auto n_samples = test_for_normality(data_sheet,options,As,n_sig);

	std::cout << data_sheet.date << ", " << std::setw(3) << data_sheet.no << "   (" << n_samples << " samples)\n";
	/*std::cout << std::setprecision(3);
	std::cout << "A²* = {";
	output_separated(std::cout,As.begin(),As.end(),", ");
	std::cout << "}\n";*/

	auto n_50 = n_sig.at(0);
	auto n_25 = n_sig.at(1);
	auto n_10 = n_sig.at(3);
	auto n_5  = n_sig.at(4);
	auto n_1  = n_sig.at(5);

	bool space = false;
	std::cout << "Level   :     50%      25%      10%       5%       1%\n";
	std::cout << "Rejected: ";
	std::cout << std::setw(7) << n_of_m(n_50,N,space) << "  " <<
				 std::setw(7) << n_of_m(n_25,N,space) << "  " <<
				 std::setw(7) << n_of_m(n_10,N,space) << "  " <<
				 std::setw(7) << n_of_m(n_5 ,N,space) << "  " <<
				 std::setw(7) << n_of_m(n_1 ,N,space) << "\n";


	std::cout << "\n";
}





void execute_filename(const std::string& filename,const DataSheet& data_sheet,const Options& options)
{
	if (options.stats) {
		write_data_sheet_stats(std::cout,data_sheet);
		std::cout << "\n";
	}
	else {
		std::cout << filename << "\n";
	}
}



void execute_header(const DataSheet& data_sheet,const Options& )
{
	write_header(std::cout,data_sheet);
}



void execute_raw(const DataSheet& data_sheet,const Options& )
{
	write_raw(std::cout,data_sheet);
}



void execute_reduce_raw(const DataSheet& data_sheet,const Options& options)
{
	std::vector<ReducedTurn> reduced_turns;

	selected_and_transformed_turns(data_sheet,options,[&](int i,const DataSheet::Turn& turn,const std::array<float,17>& distances){
		auto displacements = reduce_from_drift_and_offset(distances);

		if (options.single)
			reduce_to_single_period(displacements);

		reduced_turns.push_back({i,turn,std::move(displacements)});
	});


	write_reduce_raw(std::cout,reduced_turns);
}



std::map<double,std::complex<double>> 
raw_spectrum(const DataSheet& data_sheet,const Options& options)
{
	const int sample_rate = 16; // per turn
	
	std::vector<double> frequencies;
	int n=0;
	selected_and_transformed_turns(data_sheet,options,[&n](int,const DataSheet::Turn&,const std::array<float,17>& ){		
		n++;
	});
	// calculate valid frequencies
	for (int i=1;i<n*sample_rate/2;i++) { // skip 0
		frequencies.push_back(i/double(n));
	}
	
	// This is not efficient, but do it anyway!
	DFTGoertzel dft(frequencies,sample_rate);	
	float offset = 0;
	selected_and_transformed_turns(data_sheet,options,[&dft,&offset](int,const DataSheet::Turn&,const std::array<float,17>& distances){				
		auto ds = distances;
		for (auto& d : ds)
			d -= distances[0] - offset; // de-adjust, remove jumps
		dft.analyze(ds.begin(),ds.end()-1);
		offset += distances[16]-distances[0];
	});
	
	return dft.result();
}



void execute_raw_spectrum(const DataSheet& data_sheet,const Options& options)
{
	auto result = raw_spectrum(data_sheet,options);
	
	write_spectrum(std::cout,result);
}



std::map<double,std::complex<double>> 
spectrum(const DataSheet& data_sheet,const Options& options)
{
	const int sample_rate = 16; // per turn
	
	std::vector<double> frequencies;
	int n=0;
	selected_and_transformed_turns(data_sheet,options,[&n](int,const DataSheet::Turn&,const std::array<float,17>& ){		
		n++;
	});
	// calculate valid frequencies
	for (int i=1;i<n*sample_rate/2;i++) { // skip 0
		frequencies.push_back(i/double(n));
	}
	
	// This is not efficient, but do it anyway!
	DFTGoertzel dft(frequencies,sample_rate);	
	selected_and_transformed_turns(data_sheet,options,[&dft,&options](int,const DataSheet::Turn&,const std::array<float,17>& distances){
		auto displacements = reduce_from_drift_and_offset(distances);

		if (options.single)
			reduce_to_single_period(displacements);
		// different offsets between turns are ignored, they are small
		dft.analyze(displacements.begin(),displacements.end()-1);				
	});
	
	return dft.result();
}



void execute_spectrum(const DataSheet& data_sheet,const Options& options)
{
	
	auto result = spectrum(data_sheet,options);
	
	write_spectrum(std::cout,result);
}



void execute(Action action,const Options& options,const DataSheet& data_sheet,const std::string& filename)
{

	switch(action) {
	case Action::Filename:
		execute_filename(filename,data_sheet,options);
		break;		
	case Action::Header:
		execute_header(data_sheet,options);
		break;
	case Action::Raw:
		execute_raw(data_sheet,options);
		break;
	case Action::RawReduced:
		execute_reduce_raw(data_sheet,options);
		break;
	case Action::Reduce:
		execute_reduce(data_sheet,options);
		break;		
	case Action::Test:
		execute_test(data_sheet,options);
		break;
	case Action::Spectrum:
		execute_spectrum(data_sheet,options);
		break;
	case Action::RawSpectrum:
		execute_raw_spectrum(data_sheet,options);
		break;
	default:
		throw std::runtime_error("unknown action");
	}
}





struct MetaData
{
	bool normal = false;
	bool inverted = false;
	bool desk_in_sw = false;
	bool desk_in_nw = false;
	bool day = false;
	bool night = false;
	bool epoch_apr = false;
	bool epoch_aug = false;
	bool epoch_sep = false;
	bool epoch_feb = false;
};



void collect_meta_data(const DataSheet& data_sheet,MetaData& meta_data)
{
	meta_data.normal = meta_data.normal || !data_sheet.inverted;
	meta_data.inverted = meta_data.inverted || data_sheet.inverted;
	meta_data.desk_in_nw = meta_data.desk_in_nw || !data_sheet.desk_in_sw;
	meta_data.desk_in_sw = meta_data.desk_in_sw || data_sheet.desk_in_sw;
	auto td = time_of_day(data_sheet);
	meta_data.day = meta_data.day || td!=TimeOfDay::Night;
	meta_data.night = meta_data.night || td==TimeOfDay::Night;
	switch(epoch(data_sheet)) {
	case Epoch::Apr:
		meta_data.epoch_apr = true;
		break;
	case Epoch::Aug:
		meta_data.epoch_aug = true;
		break;
	case Epoch::Sep:
		meta_data.epoch_sep = true;
		break;
	case Epoch::Feb:
		meta_data.epoch_feb = true;
		break;
	default:
		throw std::runtime_error("unknown epoch in collect meta data");
	}
}



void validate_meta_data(const MetaData& meta_data)
{
	int epochs=0;
	if (meta_data.epoch_apr)
		epochs++;
	if (meta_data.epoch_aug)
		epochs++;
	if (meta_data.epoch_sep)
		epochs++;
	if (meta_data.epoch_feb)
		epochs++;

	if (epochs > 1) {
		std::cerr << "WARNING: Aggregating data sheets from different epochs\n";
	}
	if (meta_data.normal && meta_data.inverted) {
		std::cerr << "WARNING: Aggregating normal and inverted data sheets\n";
	}

	if (meta_data.desk_in_nw && meta_data.desk_in_sw) {
		std::cerr << "WARNING: Aggregating NW and SW data sheets\n";
	}

	if (meta_data.day && meta_data.night) {
		std::cerr << "WARNING: Aggregating day and night data sheets\n";
	}
}



DisplacementData
aggregate_reduce(const std::vector<std::reference_wrapper<const DataSheet>>& data_sheets,const Options& options)
{	
	DisplacementData aggregated_displacements {};	
	
	MetaData meta_data;
	for (const DataSheet& data_sheet : data_sheets) {
		collect_meta_data(data_sheet,meta_data);		
		auto displacements = process_data_sheet(data_sheet,options);

		add_array(aggregated_displacements.data,displacements.data);
		add_sqr_array(aggregated_displacements.uncertainties,displacements.uncertainties);
		add_array(aggregated_displacements.theory,displacements.theory);
		add_array(aggregated_displacements.model,displacements.model);
										
		aggregated_displacements.z1.m += displacements.z1.m;		
		aggregated_displacements.z2.m += displacements.z2.m;
		aggregated_displacements.z1.u += std::complex<double> {sqr(displacements.z1.u.real()), sqr(displacements.z1.u.imag())};
		aggregated_displacements.z2.u += std::complex<double> {sqr(displacements.z2.u.real()), sqr(displacements.z2.u.imag())};
	}
	validate_meta_data(meta_data);
	
	for (auto& distance : aggregated_displacements.data)
		distance /= data_sheets.size();
	for (auto& u : aggregated_displacements.uncertainties)
		u = std::sqrt(u)/data_sheets.size();
	for (auto& distance : aggregated_displacements.theory)
		distance /= data_sheets.size();
	for (auto& distance : aggregated_displacements.model)
		distance /= data_sheets.size();
	
		
	aggregated_displacements.z1.m /= data_sheets.size();
	aggregated_displacements.z2.m /= data_sheets.size();
	auto uR1 = std::sqrt(aggregated_displacements.z1.u.real())/data_sheets.size();
	auto uI1 = std::sqrt(aggregated_displacements.z1.u.imag())/data_sheets.size();
	aggregated_displacements.z1.u = {uR1,uI1};		
	auto uR2 = std::sqrt(aggregated_displacements.z2.u.real())/data_sheets.size();
	auto uI2 = std::sqrt(aggregated_displacements.z2.u.imag())/data_sheets.size();
	aggregated_displacements.z2.u = {uR2,uI2};		

	return aggregated_displacements;
}



void execute_aggregate_mean_reduce(const std::vector<DataSheet>& data_sheets,const Options& options)
{
	std::vector<std::reference_wrapper<const DataSheet>> data_sheets_view;
	for (const auto& data_sheet : data_sheets) {
		data_sheets_view.push_back(std::ref(data_sheet));
	}

	auto displacement_data = aggregate_reduce(data_sheets_view,options);
	write_aggregated_data(std::cout,displacement_data,options);
}





/**
 * \~german
 * Amplitude schätzen
 *
 * Versucht die Amplitude zu schätzen.
 *
 * \~english
 * Try to estimate the amplitude.
 *
 * \~
 * @param data
 * @return
 */
double estimate_amplitude(const std::array<double,17>& data)
{	
	auto mm = std::minmax_element(data.begin(),data.end());
	return (std::abs(*mm.first) + std::abs(*mm.second))*0.5;
}






/**
 *
 * @param sidereal_time in h
 * @return 
 */
int sidereal_time_to_key(double sidereal_time)
{
	return std::floor(sidereal_time/Sidereal_Aggregation_Bin_Width);
}

double key_to_mean_sidereal_time(int key)
{
	return key*Sidereal_Aggregation_Bin_Width + Sidereal_Aggregation_Bin_Width/2;
}



void execute_aggregate_sidereal(const std::vector<DataSheet>& data_sheets,const Options& options)
{
	std::map<int,SiderealData> aggregated_sidereal_data;

	for (const DataSheet& data_sheet : data_sheets) {				
		auto displacement = process_data_sheet(data_sheet,options);

		SiderealData sidereal_data;
		sidereal_data.data_amplitude = estimate_amplitude(displacement.data);
		sidereal_data.theory_amplitude = max_abs_value(displacement.theory);
		sidereal_data.model_amplitude = max_abs_value(displacement.model);

		DataSheetStats stats = data_sheet_stats(data_sheet);
		if (stats.mean_T.has_value()) {
			sidereal_data.TD = *stats.max_TD;			
		}

		int key = sidereal_time_to_key(time_to_h(data_sheet.sidereal_time));
		auto it = aggregated_sidereal_data.find(key);
		if (it == aggregated_sidereal_data.end()) {
			aggregated_sidereal_data.insert({key, sidereal_data});
		}
		else {
			it->second.data_amplitude += sidereal_data.data_amplitude;
			it->second.theory_amplitude += sidereal_data.theory_amplitude;
			it->second.model_amplitude += sidereal_data.model_amplitude;			
			it->second.n++;

			if (stats.mean_T.has_value()) {				
				it->second.TD += sidereal_data.TD;
				it->second.nT++;
			}
		}

	}

	for (auto& sd : aggregated_sidereal_data) {
		sd.second.data_amplitude /= sd.second.n;
		sd.second.theory_amplitude /= sd.second.n;
		sd.second.model_amplitude /= sd.second.n;

		sd.second.TD /= sd.second.nT;
	}


	write_aggregated_data(std::cout,aggregated_sidereal_data,options);
}




void execute_aggregate_diff(std::vector<DataSheet>& data_sheets,const Options& options)
{
	if (data_sheets.size() < 2) {
		std::cerr << "At least two data sheets are required\n";
		throw ExitException();
	}

	std::sort(data_sheets.begin(),data_sheets.end());

	auto current_sheet = data_sheets.cbegin();
	auto last_sheet = current_sheet;	
	auto last_displs = process_data_sheet(*last_sheet,options);

	std::vector<double> rs;
	for (++current_sheet; current_sheet!=data_sheets.cend(); ++current_sheet) {				
		auto displs = process_data_sheet(*current_sheet,options);
		double R = 0;
		auto d0i = last_displs.data.begin();
		auto u0i = last_displs.uncertainties.begin();
		auto d1i = displs.data.begin();
		auto u1i = displs.uncertainties.begin();
		// ignore options like -single
		for (;d0i!=last_displs.data.end(); ++d0i,++u0i,++d1i,++u1i) {
			R += sqr((*d1i - *d0i)/std::sqrt(sqr(*u0i)+sqr(*u1i))); // squared residuals with error propagation
		}
		rs.push_back(R);
		last_sheet = current_sheet;
		last_displs = displs;
	}

	std::cout << std::setprecision(2) << std::fixed;
	std::cout << "{";
	output_separated(std::cout,rs.begin(),rs.end(),", ");
	std::cout << "}\n";

	auto mean_r = mean_value(rs.begin(),rs.end());
	std::cout << "mean R² = " << mean_r << "\n";
	std::cout << "\n";
}




struct FitModelParameters
{
	double chi;
	ErrorModelParameters params;
};

bool operator >(const FitModelParameters& a,const FitModelParameters& b)
{
	return a.chi > b.chi;
}



std::vector<FitModelParameters>
find_model_parameters(const std::vector<DataSheet>& data_sheets,const Options& options)
{
	using reversed_priority_queue = std::priority_queue<FitModelParameters,std::vector<FitModelParameters>,std::greater<FitModelParameters>>;
	std::map<double,std::map<int,reversed_priority_queue>> best_params;


	std::vector<DisplacementData> displacement_data;
	MetaData meta_data;
	for (const DataSheet& data_sheet : data_sheets) {
		collect_meta_data(data_sheet,meta_data);
		displacement_data.push_back(process_data_sheet(data_sheet,options));
	}
	validate_meta_data(meta_data);

	const bool desk_in_sw = data_sheets.front().desk_in_sw;
	const int start_azimuth = desk_in_sw ? 10 : 14; // use first, should be all the same
	for (double desk_azimuth=start_azimuth-2; desk_azimuth<=start_azimuth+2; desk_azimuth+=0.25) {
		for (double desk=0; desk<0.2; desk+=0.005) {
			for (double amplitude=0; amplitude<0.2; amplitude+=0.005) {
				ErrorModelParameters params {desk_azimuth,desk,amplitude};

				double chi = 0;
				for_each_of(displacement_data,data_sheets, [&](const DisplacementData& displacements,const DataSheet& data_sheet) {
					auto model = systematic_error_displacements(data_sheet.thermometers_start,data_sheet.thermometers_end,
																params,options);
					if (!options.subtract)
						add_array(model,displacements.theory);

					chi += chi_squared_test(displacements.data,displacements.uncertainties,model);
				});

				chi /= data_sheets.size();
				auto groupId = std::abs(desk-amplitude)<0.01 ? 0 : (desk<amplitude ? -1 : 1);
				best_params[desk_azimuth][groupId].push({chi,params});
			}
		}
	}

	std::vector<FitModelParameters> ordered_best_params;
	for (auto& azimuth_groups : best_params) {
		for (auto& group_params : azimuth_groups.second) {
			int n=0;
			while (!group_params.second.empty() && n<1) {
				ordered_best_params.push_back(group_params.second.top());
				group_params.second.pop();
				n++;
			}
		}
	}

	return ordered_best_params;
}



void write_aggregated_data(std::ostream& os,const std::vector<FitModelParameters>& params)
{

	for (const FitModelParameters& param : params) {
		os << std::setprecision(2) << std::fixed;
		os << "χ² = " << std::setw(5) << param.chi << ", params: {";
		os << std::setprecision(3) << std::fixed;
		output_args_separated(os,", ",param.params.desk_azimuth,
							  param.params.desk,param.params.temperature);
		os << "}\n";
	}
	os << "\n";
}



void execute_aggregate_params(const std::vector<DataSheet>& data_sheets,const Options& options)
{
	auto found_params = find_model_parameters(data_sheets,options);

	write_aggregated_data(std::cout,found_params);
}




void execute_aggregate_model_chi(const std::vector<DataSheet>& data_sheets,const Options& options)
{
	if (!options.model) {
		std::cerr << "Model must be activated\n";
		throw ExitException();
	}

	const int N = azimuths(options);

	double chi = 0.0;
	//MetaData meta_data;
	for (const auto& data_sheet : data_sheets) {
		//collect_meta_data(data_sheet,meta_data);
		auto displacement_data = process_data_sheet(data_sheet,options);
		chi += chi_squared_test(displacement_data.data, displacement_data.uncertainties, displacement_data.model,N);
	}
	//validate_meta_data(meta_data);


	auto n = data_sheets.size()*N;
	// Estimate number of Parameters
	if (options.epoch_params) {
		n -= 3*data_sheets.size()/40; // 1 epoch ~80, with distinction sw/nw ~40
	}
	else {
		// 3 model params, group size ~4
		n -= 3*data_sheets.size()/4;
	}
	auto p = 1.-chi_square_cdf(chi,n);

	write_chi_squared_stats(std::cout,chi,n,p);
	std::cout << "The degrees of freedom f are just estimated!\n";
}



void execute_aggregate_list(std::vector<DataSheet>& data_sheets,const Options& options)
{
	std::sort(data_sheets.begin(),data_sheets.end());

	write_list(std::cout,data_sheets,options);
}



/**
 * \~german Ablehnungsquote
 * \~english
 * \~
 * @param r number of rejected test
 * @param n number of tests
 * @return 
 */
ConfidenceInterval rejection_quota(int r,int n)
{	
	return confidence_interval(r,n);	
}


std::string to_string(const ConfidenceInterval& q)
{
	std::ostringstream oss;
	oss << std::setprecision(1) << std::fixed;
	oss << q.p*100 << " ± " << q.u*100 << " %";
	return oss.str();
}



void execute_aggregate_test(const std::vector<DataSheet>& data_sheets,const Options& options)
{
	const auto N = azimuths(options);
	size_t total_samples = 0;
	std::vector<int> total_sig(6);

	std::vector<double> As;
	std::vector<int> n_sig;
	for (auto& data_sheet : data_sheets) {
		As.clear();
		n_sig.clear();
		total_samples += test_for_normality(data_sheet,options, As,n_sig);
		std::transform(n_sig.begin(),n_sig.end(),total_sig.begin(),total_sig.begin(),std::plus<int>());
	}

	auto mean_samples = double(total_samples)/data_sheets.size();

	std::cout << std::setprecision(1) << std::fixed;
	std::cout << "Mean number of samples per azimuth: " << mean_samples << "\n";

	auto n_50 = total_sig.at(0);
	auto n_25 = total_sig.at(1);
	auto n_10 = total_sig.at(3);
	auto n_5  = total_sig.at(4);
	auto n_1  = total_sig.at(5);

	int n = N*data_sheets.size();
	std::cout << " Level    Rejected      Quota\n";
	std::cout << "  50% "
			  << std::setw(12) << n_of_m(n_50,n) << "  "
			  << std::setw(14) << to_string(rejection_quota(n_50,n)) << "\n";
	std::cout << "  25% "
			  << std::setw(12) << n_of_m(n_25,n) << "  "
			  << std::setw(14) << to_string(rejection_quota(n_25,n)) << "\n";
	std::cout << "  10% "
			  << std::setw(12) << n_of_m(n_10,n) << "  "
			  << std::setw(14) << to_string(rejection_quota(n_10,n)) << "\n";
	std::cout << "   5% "
			  << std::setw(12) << n_of_m(n_5,n) << "  "
			  << std::setw(14) << to_string(rejection_quota(n_5,n)) << "\n";
	std::cout << "   1% "
			  << std::setw(12) << n_of_m(n_1,n) << "  "
			  << std::setw(14) << to_string(rejection_quota(n_1,n)) << "\n";


}



std::string remove_comments(const std::string& line)
{
	auto tokens = split(line,'#');
	return tokens.at(0);
}

void output_expression(std::ostream& os,const SignalExtractionExpression& expr)
{		
	os << epoch(expr.epoch.from) << " ";
		
	output_separated(os,expr.left.begin(),expr.left.end(),",",[](IntegerInterval interval){
		std::ostringstream oss;
		oss << "[";
		oss << interval.from << "," << interval.to;
		oss << "]";
		return oss.str();
	});
	
	os << " - ";
	
	output_separated(os,expr.right.begin(),expr.right.end(),",",[](IntegerInterval interval){
		std::ostringstream oss;
		oss << "[";
		oss << interval.from << "," << interval.to;
		oss << "]";
		return oss.str();
	});

}

void output_expressions(std::ostream& os,const std::vector<SignalExtractionExpression>& expressions)
{
	for (auto& expr : expressions) {
		output_expression(os,expr);
		os << "\n";
	}
}

bool are_overlapping(const IntegerInterval& a, const IntegerInterval& b)
{
	return (a.from>=b.from && a.from<=b.to) || (a.to>=b.from && a.to<=b.to);
}

bool validate_fit_expression(const std::vector<IntegerInterval>& sub_expr)
{
	if (sub_expr.size()<2)
		return true;

	for (auto i1 = sub_expr.begin();i1+1!=sub_expr.end();++i1) {
		for (auto i2 = i1+1;i2!=sub_expr.end();++i2) {
			if (are_overlapping(*i1,*i2)) {
				return false;
			}
		}
	}
	return true;
}

bool validate_fit_expression(const SignalExtractionExpression& expr)
{
	for (auto i1 = expr.left.begin();i1!=expr.left.end();++i1) {
		for (auto i2 = expr.right.begin();i2!=expr.right.end();++i2) {
			if (are_overlapping(*i1,*i2)) {
				return false;
			}
		}
	}
	return true;
}

bool validate_fit_expression(std::ostream& os,const SignalExtractionExpression& expr)
{
	if (!validate_fit_expression(expr.left)) {
		os << "\nWARNING: Overlapping intervals or duplicate numbers on left side\n";
		return false;
	}
	if (!validate_fit_expression(expr.right)) {
		os << "\nWARNING: Overlapping intervals or duplicate numbers on right side\n";
		return false;
	}
	if (!validate_fit_expression(expr)) {
		os << "\nWARNING: Overlapping intervals or duplicate numbers\n";
		return false;
	}
	
	return true;
}


Filter create_filter(const IntegerInterval& months,const std::vector<IntegerInterval>& sheet_numbers)
{
	Filter filter;
	filter.month.emplace_back(months.from, months.to);
	for (const auto& nos : sheet_numbers) {
		filter.no.emplace_back(nos.from, nos.to);
	}
	return filter;
}



struct ExtractedSignal
{
	std::vector<std::reference_wrapper<const DataSheet>> left_data_sheets;
	std::vector<std::reference_wrapper<const DataSheet>> right_data_sheets;
	std::array<double,17> data;
	std::array<double,17> uncertainties;	
	SignalParameters<Estimate<double>> k2;	
};



SignalParameters<Estimate<double>> 
extract_signal_parameters(const std::array<double,17>& data,const std::array<double,17>& uncertainties,
						  const Options& options)
{
	SignalParameters<Estimate<double>> params;
	
	auto xsine = fit_sine(data,uncertainties,options);	
	if (!xsine.valid) {
		//std::cerr << "Failed to converge at sine fitting\n";					
		//throw ExitException();
	}
	
	// the uncertainties are probably not correct, but not fully wrong
	params.A.m = xsine.x[0];
	params.A.u = xsine.u[0];
	params.phi.m = xsine.x[1];
	params.phi.u = xsine.u[1];
	
	return params;
}



SignalParameters<Estimate<double>> 
extract_signal_parameters(const Estimate<std::complex<double>>& lz2,const Estimate<std::complex<double>>& rz2)
{
	SignalParameters<Estimate<double>> params;
	
	auto z = propagate_sub(lz2, rz2); // interference
			
	const double h = 1e-6;
	params.A = propagate_uncertainties(Real(z),Imag(z),[](double R,double I){
		return std::abs(std::complex<double>{R,I});	
	},h);
	
	params.phi = propagate_uncertainties(Real(z),Imag(z),[](double R,double I){
		return std::arg(std::complex<double>{R,I});	
	},h);					
		
	return params;
}



std::vector<ExtractedSignal> extract_signals(const std::vector<SignalExtractionExpression>& expressions,
											 const std::vector<DataSheet>& data_sheets,const Options& options)
{
	std::vector<ExtractedSignal> extracted_signals;

	Options expr_options = options;
	expr_options.aggregation_method = Options::AggregationMethod::Mean;
	expr_options.model = false;

	Options left_expr_options = expr_options;
	left_expr_options.subtract_data = true;

	int i=0;
	for (const auto& expr : expressions) {
		i++;
		if (options.disabled_signals.find(i)!=options.disabled_signals.end())
			continue;
		ExtractedSignal extracted_signal;			

		// aggregate right side of expression
		Filter right_filter = create_filter(expr.epoch,expr.right);
		for (const auto& data_sheet : data_sheets) {
			if (selected(data_sheet,expr_options,right_filter)) {
				extracted_signal.right_data_sheets.push_back(std::ref(data_sheet));
			}			
		}

		auto right_agg_displs = aggregate_reduce(extracted_signal.right_data_sheets,expr_options);

		// set result of right aggregation as data for subtraction
		left_expr_options.data1.clear();
		left_expr_options.data2.clear();
		left_expr_options.data3.clear();
		std::copy(right_agg_displs.data.begin(),right_agg_displs.data.end(),
				  std::back_inserter(left_expr_options.data1));
		std::copy(right_agg_displs.uncertainties.begin(),right_agg_displs.uncertainties.end(),
				  std::back_inserter(left_expr_options.data2));
		std::copy(right_agg_displs.theory.begin(),right_agg_displs.theory.end(),
				  std::back_inserter(left_expr_options.data3));

		// aggregate left side of expression
		Filter left_filter = create_filter(expr.epoch,expr.left);
		for (const auto& data_sheet : data_sheets) {
			if (selected(data_sheet,left_expr_options,left_filter)) {
				extracted_signal.left_data_sheets.push_back(std::ref(data_sheet));
			}			
		}

		if (extracted_signal.left_data_sheets.empty() || extracted_signal.right_data_sheets.empty()) {
			std::cerr << "skipping signal: ";
			output_expression(std::cerr,expr);
			std::cerr << "\n";
			continue;
		}

		auto left_agg_displs = aggregate_reduce(extracted_signal.left_data_sheets,left_expr_options);

		// add the resulting signal
		extracted_signal.data = left_agg_displs.data;
		extracted_signal.uncertainties = left_agg_displs.uncertainties;
		if (options.fit_sine) {
			if (options.reduction_method == Options::ReductionMethod::DFT) {
				extracted_signal.k2 = extract_signal_parameters(left_agg_displs.z2,right_agg_displs.z2);		
			}
			else {
				extracted_signal.k2 = extract_signal_parameters(left_agg_displs.data,left_agg_displs.uncertainties,options);		
			}
		}
				
		extracted_signals.push_back(extracted_signal);
	}

	return extracted_signals;
}



std::array<double,17> theory_signal(const Theory& theory,const TheoryParameters& params,
									const ExtractedSignal& extracted_signal,const Options& options)
{
	// left
	std::array<double,17> left_displs {};
	for (const DataSheet& data_sheet : extracted_signal.left_data_sheets) {		
		auto displacements = fringe_displacements(theory,params,data_sheet,options);
		add_array(left_displs,displacements);
	}
	for (auto& displ : left_displs)
		displ /= extracted_signal.left_data_sheets.size(); // mean values

	// right
	std::array<double,17> right_displs {};
	for (const DataSheet& data_sheet : extracted_signal.right_data_sheets) {		
		auto displacements = fringe_displacements(theory,params,data_sheet,options);
		add_array(right_displs,displacements);
	}
	for (auto& displ : right_displs)
		displ /= extracted_signal.right_data_sheets.size(); // mean values

	// subtract
	sub_array(left_displs,right_displs);
	return left_displs;
}




struct SignalStats {
	std::vector<double> chis;
	std::vector<double> residuals;		
};



const std::array<double,17> theory_displs_u {1,1,1,1, 1,1,1,1, 1,1,1,1, 1,1,1,1, 1};



double chi_squared_sum(const Theory& theory, const TheoryParameters& params,
					   const std::vector<ExtractedSignal>& extracted_signals,const Options& options,
					   bool collect_stats,SignalStats& signal_stats)
{
	double chi_sum = 0;
	const auto N = azimuths(options);
		
	for (const auto& extracted_signal : extracted_signals) {
		double chi;
		auto theory_displs = theory_signal(theory,params,extracted_signal,options);		
		if (options.fit_amplitude) {
			double R;
			if (options.fit_sine) {		
				if (options.reduction_method == Options::ReductionMethod::DFT) {
					auto tz = DFT_analyze(2,theory_displs.begin(),theory_displs.end()-1);
					auto A = std::abs(tz);
					
					R = (extracted_signal.k2.A.m-A)/extracted_signal.k2.A.u;
				}
				else {
					auto tsine = fit_sine(theory_displs,theory_displs_u,options);
					if (!tsine.valid) { // seems to work anyway
						//std::cerr << "Failed to converge at sine fitting\n";					
						//throw ExitException();
					}
								
					R = (extracted_signal.k2.A.m-tsine.x[0])/extracted_signal.k2.A.u;
				}
				chi = sqr(R);											
			}
			else {
				auto signal_amp = estimate_amplitude(extracted_signal.data);
				auto amp_u = mean_value(extracted_signal.uncertainties)*0.5; // estimate; Empiric "proof" in fit_sine test case.
				auto theory_amp = *std::max_element(theory_displs.begin(),theory_displs.end());
							
				R = (signal_amp-theory_amp)/amp_u;
				chi = sqr(R);					
			}
			
			if (collect_stats) {				
				signal_stats.chis.push_back(chi);											
				signal_stats.residuals.push_back(R);
			}
		}
		else {		
			if (options.fit_sine) {				
				double Rp,Ra;
				if (options.reduction_method == Options::ReductionMethod::DFT) {
					auto tz = DFT_analyze(2,theory_displs.begin(),theory_displs.end()-1);
					auto A = std::abs(tz);
					auto phi = std::arg(tz);
					
					auto pd = phase_difference(phi,extracted_signal.k2.phi.m,AETHER_2PI); // cosine phase
					Rp = pd/extracted_signal.k2.phi.u;
					Ra = (extracted_signal.k2.A.m-A)/extracted_signal.k2.A.u;						
				}
				else {					
					auto tsine = fit_sine(theory_displs,theory_displs_u,options);
					if (!tsine.valid) {
						//std::cerr << "Failed to converge at sine fitting\n";
						//throw ExitException();
					}	
															
					auto pd = phase_difference(tsine.x[1],extracted_signal.k2.phi.m,AETHER_2PI); // sine phase
					Rp = pd/extracted_signal.k2.phi.u;
					Ra = (extracted_signal.k2.A.m-tsine.x[0])/extracted_signal.k2.A.u;													
				}
				
				chi = sqr(Rp) + sqr(Ra); // phase + amplitude
								
				if (collect_stats) {						
					signal_stats.residuals.push_back(Rp);
					signal_stats.residuals.push_back(Ra);
				}
			}
			else {
				chi = chi_squared_test(extracted_signal.data,extracted_signal.uncertainties,theory_displs,N);
				if (collect_stats) {					
					auto rs = residuals(extracted_signal.data,extracted_signal.uncertainties,theory_displs,N);					
					for (auto ri = rs.begin(); ri!=rs.begin()+N; ++ri)
						signal_stats.residuals.push_back(*ri);					
				}
			}
			
			if (collect_stats) {				
				signal_stats.chis.push_back(chi);			
			}
		}	
		chi_sum += chi;
	}	

	return chi_sum;
}



double chi_squared_sum(const Theory& theory, const TheoryParameters& params,
					   const std::vector<ExtractedSignal>& extracted_signals,const Options& options) {
	SignalStats stats;
	return chi_squared_sum(theory,params,extracted_signals,options,false,stats);
}



/**
 * \~german
 * Unsicherheit der Parameter am Minimum. 
 * Ohne Beachtung von Korrelationen.
 * 
 * \~english
 * Uncertainty of the parameters at the minimum.
 * Without taking correlations into account.
 * 
 * \~
 * @param min_params parameters at the minimum
 * @param min_chi chi squared value at the minimum
 * @param dchi2 delta chi squared value to search for
 * @param min minimum parameter value (usefull for periodic parameters)
 * @param max maximum parameter value (usefull for periodic parameters)
 * @param f chi squared function object
 * @param h step
 * @param temp_params used internaly
 * @param selected_temp_param used internaly
 * @return uncertainty
 */
template<typename X,typename T=typename X::value_type,typename F>
double fit_parameter_uncertainty(const X& min_params,double min_chi,double dchi2,double min,double max,
								 F f,double h, X& temp_params,T& selected_temp_param)
{		
	// TODO This is too simple, hypercontour should be scanned.
	// Errors are correct only for small covariances/correlations.
	temp_params = min_params;
	double min_param = selected_temp_param;

	double chi;
	double c = h; // can be negative
	double last_param;
	do {
		last_param = selected_temp_param;
		selected_temp_param = min_param+c;					
		chi = f(temp_params);
		if (c<min || c>max) { // prevent loop for periodic params
			c = clamp(selected_temp_param,min,max);
			selected_temp_param = min_param+c;
			break;
		}
		c+=c;
	} 
	while(chi-min_chi<dchi2 && !std::isinf(c));

	if (chi-min_chi>dchi2) {
		// approach chi = dchi2, which lies between last_param and selected_temp_param
		bisection(last_param,selected_temp_param,[&](double ,double x,double ) {
			selected_temp_param = x;
			auto chi = f(temp_params);
			auto d = chi - min_chi;
			if (std::abs(d-dchi2)<0.000001)
				return 0.0; // cancel, close enough
			return dchi2-d; // choose direction and continue
		});
	}

	return std::abs(selected_temp_param-min_param);
}



double fit_parameter_uncertainty(const Theory& theory,const TheoryParameters& min_params,double min_chi,double h,double min,double max,
								 const std::vector<ExtractedSignal>& extracted_signals,const Options& options,
								 TheoryParameters& temp_params,double& selected_temp_param)
{
	auto f = [&theory,&extracted_signals,&options](const TheoryParameters& params){
		return chi_squared_sum(theory, params,extracted_signals,options);
	};
	return fit_parameter_uncertainty(min_params,min_chi,options.delta_chi_squared,min,max,f,h,temp_params,selected_temp_param);
}



TheoryParameters
fit_parameters_uncertainties(const TheoryParameters& min_params,double min_chi,
							 const std::vector<ExtractedSignal>& extracted_signals,const Options& options)
{
	TheoryParameters u;
	TheoryParameters temp_params;
	auto theory = create_theory(options.theory);
	const auto inf = std::numeric_limits<double>::infinity();

	double neg_v = fit_parameter_uncertainty(*theory,min_params,min_chi,-1,-inf,inf,extracted_signals,options,
											 temp_params,temp_params.v);
	double pos_v = fit_parameter_uncertainty(*theory,min_params,min_chi,1,-inf,inf,extracted_signals,options,
											 temp_params,temp_params.v);
	u.v = std::max(neg_v,pos_v);


	double neg_ra = fit_parameter_uncertainty(*theory,min_params,min_chi,-0.000001,h_to_rad(-12),h_to_rad(12),
											  extracted_signals,options, temp_params,temp_params.a);
	double pos_ra = fit_parameter_uncertainty(*theory,min_params,min_chi,0.000001,h_to_rad(-12),h_to_rad(12),
											  extracted_signals,options, temp_params,temp_params.a);
	u.a = std::max(neg_ra,pos_ra);


	double neg_de = fit_parameter_uncertainty(*theory,min_params,min_chi,-0.000001,rad(-90),rad(90),
											  extracted_signals,options, temp_params,temp_params.d);
	double pos_de = fit_parameter_uncertainty(*theory,min_params,min_chi,0.000001,rad(-90),rad(90),
											  extracted_signals,options, temp_params,temp_params.d);
	u.d = std::max(neg_de,pos_de);


	return u;
}



/**
 * \~german
 * Fortschrittsanzeige in Prozent
 *
 * Ausgabe und Aktualisierung auf der Konsole an der aktuellen
 * Zeigerposition.
 *
 * \~english
 * Display progress in percent
 *
 * Output and update on the console at the current cursor position.
 */
class ProgressInPercent
{
	int last_percent = -1;
	std::ostream& os;
public:
	ProgressInPercent(std::ostream& os)
		:os {os}{
		os << store_cursor;
	}

	void operator() (int i,int n) {
		int percent = int((100.0/n)*i);
		if (percent!=last_percent) {
			os << restore_cursor << percent << "%" << std::flush;
			last_percent = percent;
		}
	}

};



class EstimateProgress
{		
	int max_iterations;
	double log_prec;
	int last_j = 0;
	ProgressInPercent progress;

public:
	EstimateProgress(std::ostream& os,double precision,int max_iterations)
		:max_iterations{max_iterations},progress {os}
	{
		log_prec = std::log10(precision);
	}

	void operator() (int i,double y,double y1) {		
		auto log_y = std::log10(std::abs(1-y1/y));
		int j = max_iterations/log_prec*log_y;
		j = clamp(j,last_j,max_iterations);

		progress(std::max(i,j),max_iterations);
		last_j = j;
	}
};



std::ostream& fit_ostream(const Options& options)
{
	return (options.contour || options.residuals) ? std::cerr : std::cout;
}



const double VScaleFactor = 100000.0;

MinimizerResult
fit_theory_grad(const TheoryParameters& params,const std::vector<ExtractedSignal>& extracted_signals,
				const Options& options)
{
	MinimizerResult result;
	std::ostream& os = fit_ostream(options);

	auto theory = create_theory(options.theory);
	const double sf = VScaleFactor; // scale factor, because h is a small value choosen for ra, de
	
	double h = 1e-5;
	int n = 100;
	double precision = 1e-8;

	
	std::array<double,3> x0 {params.v/sf,params.a,params.d};
	
	auto func = [sf,&extracted_signals,&theory,&options](const std::array<double,3>& x) {
		return chi_squared_sum(*theory,{x[0]*sf,x[1],x[2]},extracted_signals,options);
	};
	
	// Method A gets consistently near the minimum, but sometimes failes to converge
	auto local_min = minimize_locally_a(x0,func,h,n,precision,EstimateProgress(os,precision,n));

	// If method A failed to converge, try method B
	if (local_min.i==n) {
		h *= 0.1;
		os << "\nMaximum iterations reached, trying alternative..." << std::flush;
		auto lm = minimize_locally_b(local_min.x,func,h,n,precision,EstimateProgress(os,precision,n));
		if (lm.y <= local_min.y) {
			local_min = lm;			
		}
	}
	
	local_min.x[0]*=sf;

	TheoryParameters min_params {local_min.x[0],local_min.x[1],local_min.x[2]};
	auto u = fit_parameters_uncertainties(min_params,local_min.y,extracted_signals,options);
	
	result.valid = local_min.i<n;
	result.x = {local_min.x[0], local_min.x[1], local_min.x[2]};
	result.u = {u.v, u.a, u.d};
	result.y = local_min.y;
	
	
	return result;
}



#ifdef AETHER_MINUIT

class FitTheoryContext : public MinimizeContext
{
	const Theory& theory;
	const std::vector<ExtractedSignal>& extracted_signals;
	const Options& options;
	double _scale_factor = 1.0;
	double _delta_chi_squared = 1.0;
public:
	FitTheoryContext(const Theory& theory,const std::vector<ExtractedSignal>& extracted_signals,
					 const Options& options)
		:theory {theory},extracted_signals{extracted_signals},options{options} {
	}

	double operator() (const std::vector<double>& params) override
	{
		return chi_squared_sum(theory,{params.at(0)*_scale_factor,params.at(1),params.at(2)},extracted_signals,options);
	}

	double scale_factor()
	{
		return _scale_factor;
	}

	void scale_factor(double sf)
	{
		_scale_factor = sf;
	}

	double delta_chi_squared() const override {
		return _delta_chi_squared;
	}
		
	void delta_chi_squared(double delta_chi_squared)
	{
		_delta_chi_squared = delta_chi_squared;
	}
};



MinimizerResult
fit_theory_Minuit2(const TheoryParameters& params,const std::vector<ExtractedSignal>& extracted_signals,
				   const Options& options)
{
	std::ostream& os = fit_ostream(options);
	auto theory = create_theory(options.theory);
	const double sf = VScaleFactor; // scale factor, because v is big compared to ra, de
	std::vector<double> x0;
	x0.push_back(params.v/sf);
	x0.push_back(params.a);
	x0.push_back(params.d);

	FitTheoryContext context(*theory,extracted_signals,options);
	context.scale_factor(sf);
	context.delta_chi_squared(options.delta_chi_squared);
	auto result = minimize_locally_Minuit2(x0,context,os);
	if (result.x.size()>0)
		result.x.at(0)*=sf;
	if (result.u.size()>0)
		result.u.at(0)*=sf;
	return result;
}

#endif


void write_coordinates(std::ostream& os,const MinimizerResult& result)
{
	if (locale_german) {
		output_args_separated(os,"\t",quote("Signal 𝛼"), quote("Signal 𝛿"), quote(""), quote("Signal"), "\n");
	}
	else {
		output_args_separated(os,"\t",quote("Signal 𝛼"), quote("Signal 𝛿"), quote(""), quote("Signal"), "\n");
	}
	output_args_separated(os,"\t",rad_to_h(result.x[1]), deg(result.x[2]), result.x[0], NAN, "\n");
	
	write_gnuplot_data_set_separator(os);
	
	if (locale_german) {
		output_args_separated(os,"\t",quote("KHS-Dipol 𝛼"),quote("KHS-Dipol 𝛿"),quote(""),quote("KHS-Dipol"),"\n");
	}
	else {
		output_args_separated(os,"\t",quote("CMB Dipole 𝛼"),quote("CMB Dipole 𝛿"),quote(""),quote("CMB Dipole"),"\n");
	}
	output_args_separated(os,"\t",rad_to_h(CMB_dipole.a), deg(CMB_dipole.d), CMB_dipole.v, NAN, "\n");

	write_gnuplot_data_set_separator(os);
}



void contour(const MinimizerResult& result,const std::vector<ExtractedSignal>& extracted_signals,
			 const Options& options,std::ostream& os)
{
	auto theory = create_theory(options.theory);

	auto v = result.x[0];
	auto ra = result.x[1];
	auto de = result.x[2];

	auto uv = result.u[0];
	auto ura = clamp(result.u[1],0.,double(AETHER_PI));
	auto ude = clamp(result.u[2],0.,double(AETHER_PI_2));

	auto min_chi = result.y;

	// sqrt() to prevent quadratic growth of number of data points!
	const int ndp = std::ceil(std::sqrt(options.delta_chi_squared));

	std::cerr << "Scanning contour for (v,α)..." << std::flush;
	ProgressInPercent progress1(std::cerr);
	for(int i=-10*ndp;i<=10*ndp;i++) {
		for(int k=-10*ndp;k<=10*ndp;k++) {
			auto dv = i/(10.*ndp)*uv;
			auto dra = k/(10.*ndp)*ura;
			auto chi2 = chi_squared_sum(*theory,{v+dv,ra+dra,de},extracted_signals,options);
			os << v+dv << "\t";
			os << ra+dra << "\t";		
			os << chi2-min_chi << "\n";
		}
		os << "\n";
		progress1(i+10*ndp,20*ndp);
	}
	std::cerr << "\n";

	os << "\n"; // becomes gnuplot dataset separator

	std::cerr << "Scanning contour for (v,δ)..." << std::flush;
	ProgressInPercent progress2(std::cerr);
	for(int i=-10*ndp;i<=10*ndp;i++) {
		for(int k=-10*ndp;k<=10*ndp;k++) {
			auto dv = i/(10.*ndp)*uv;
			auto dde = k/(10.*ndp)*ude;
			auto chi2 = chi_squared_sum(*theory,{v+dv,ra,de+dde},extracted_signals,options);
			os << v+dv << "\t";
			os << de+dde << "\t";
			os << chi2-min_chi << "\n";
		}
		os << "\n";
		progress2(i+10*ndp,20*ndp);
	}
	std::cerr << "\n";

	os << "\n"; // becomes gnuplot dataset separator

	std::cerr << "Scanning contour for (α,δ)..." << std::flush;
	ProgressInPercent progress3(std::cerr);
	for(int i=-10*ndp;i<=10*ndp;i++) {
		for(int k=-10*ndp;k<=10*ndp;k++) {
			auto dra = i/(10.*ndp)*ura;
			auto dde = k/(10.*ndp)*ude;
			auto chi2 = chi_squared_sum(*theory,{v,ra+dra,de+dde},extracted_signals,options);
			os << ra+dra << "\t";
			os << de+dde << "\t";
			os << chi2-min_chi << "\n";
		}
		os << "\n";
		progress3(i+10*ndp,20*ndp);
	}
	std::cerr << "\n";
}



MinimizerResult
fit_theory(const TheoryParameters& start_params,const std::vector<ExtractedSignal>& extracted_signals,
		   const Options& options)
{
	MinimizerResult result;

	switch(options.minimizer) {
	case Options::Minimizer::Grad:
		result = fit_theory_grad(start_params,extracted_signals,options);
		break;
#ifdef AETHER_MINUIT
	case Options::Minimizer::Minuit2:
		result = fit_theory_Minuit2(start_params,extracted_signals,options);
		break;
#endif
	default:
		throw std::runtime_error("unknown minimizer");
	}

	return result;
}



template<typename F>
TheoryParameters best_of_random_samples(const std::vector<ExtractedSignal>& extracted_signals,
										const Options& options,F progress)
{
	double min_chi = std::numeric_limits<double>::infinity();
	TheoryParameters best_params {};
	auto theory = create_theory(options.theory);

	std::mt19937 engine(create_random_seed());
	std::uniform_real_distribution<double> dist(-1.0,1.0);

	int n1 = 10*8*16;
	int n2 = options.fit_amplitude ? n1 : 0; // amplitude-only seems harder

	for (int i=1;i<=n1;i++) {
		double v = 350000 + dist(engine)*200000;
		double ra = h_to_rad(12+dist(engine)*12);
		double de = rad(dist(engine)*89);

		TheoryParameters random_params {v,ra,de};
		auto chi = chi_squared_sum(*theory,random_params,extracted_signals,options);
		if (chi < min_chi) {
			min_chi = chi;
			best_params = random_params;
		}

		progress(i,n1+n2);
	}

	// second phase around the probable minimum
	auto min_p = best_params;
	for (int i=1;i<=n2;i++) {
		double v = min_p.v + dist(engine)*50000;
		double ra = min_p.a + h_to_rad(dist(engine)*5);
		double de = min_p.d + rad(dist(engine)*15);

		TheoryParameters random_params {v,ra,de};
		auto chi = chi_squared_sum(*theory,random_params,extracted_signals,options);
		if (chi < min_chi) {
			min_chi = chi;
			best_params = random_params;
		}

		progress(n1+i,n1+n2);
	}

	return best_params;
}




int degrees_of_freedom_Miller(int n,const Options& options)
{		
	//double N = options.fit_sine ? 2 : azimuths(options);	// naive
	
	// N was determined in simulations	
	double N; // d.o.f. per signal
	
	if (options.fit_amplitude) {
		if (options.fit_sine) {			
			N = options.single ? 2.1 : 1.1;				
		}
		else {
			N = options.single ? 2.2 : 1.4;				
		}
	}
	else {
		if (options.fit_sine) {			
			N = options.single ?  4.3 : 2;	
			N *= 1.4;			
		}
		else {
			N = options.single ? 8.3 : 19;	
			N *= 1.2;
		}				
	}
	
	//return int(n*N - 3); // 3 parameters	
	return int(std::round(n*N));
}



int degrees_of_freedom_DFT(int n,const Options& options)
{
	// N was determined in simulations (if not calculated)	
	double N; // d.o.f. per signal
	
	if (options.fit_amplitude) {
		if (options.fit_sine) {
			N = 1; // amplitude
			return int(n*N-3); // naive ansatz seems to work
		}
		else {			
			N = options.single ? 3.3 : 1.4;
		}
	}
	else {
		if (options.fit_sine) {		
			N = 2; // amplitude & phase
			return int(n*N-3); // naive ansatz seems to work
		}
		else {			
			N = options.single ? 8.9 : 21;
		}		
	}
	
	return int(std::round(n*N));
}





/**
 * \~german
 * Freiheitsgrade
 * 
 * Die Anzahl der Freiheitsgrade bei der Minimierung hängt von
 * der Anzahl der Signale und den Einstellungen ab. Außerdem
 * verhält sich die Theorie anscheinend so, daß die Freiheitsgrade
 * nicht trival zu bestimmen sind. Es wurde deshalb versucht 
 * mit Simulationen Schätzwerte zu finden.
 * 
 * \~english
 * The number of degrees of freedom at minimization depends on
 * the number of signals and the settings. Moreover, the theory 
 * seems to behave in such a way that the degrees of freedom
 * are not to be determined trivially. Therefore, an attempt 
 * was made to find estimated values with simulations.
 * 
 * \~
 * @param n number of signals
 * @param options
 * @return 
 */
int degrees_of_freedom(int n,const Options& options)
{
	switch (options.reduction_method) {
	case Options::ReductionMethod::Miller:
		return degrees_of_freedom_Miller(n,options);
	case Options::ReductionMethod::DFT:
		return degrees_of_freedom_DFT(n,options);
	default:
		throw std::runtime_error("unknown reduction method");
	}
}



void write_fit_result(std::ostream& os,const MinimizerResult& result,bool comment,int n,const Options& options)
{
	const char* cc = comment ? "# " : "";

	if (result.x.size()==3 && result.u.size()==3) {
		os << cc << "v = " << result.x[0]			<< " ± " << result.u[0]           << " m/s\n";
		os << cc << "α = " << rad_to_h(result.x[1])	<< " ± " << rad_to_h(result.u[1]) << " h\n";
		os << cc << "δ = " << deg(result.x[2])		<< " ± " << deg(result.u[2])      << " °\n";
	}
	else {
		throw std::runtime_error("invalid data in result");		
	}
	os << cc << "\n";
	
	auto f = degrees_of_freedom(n,options);	
	auto p = 1.-chi_square_cdf(result.y,f);
	write_chi_squared_stats(os,result.y,f,p,cc);	
}



#ifdef AETHER_MINUIT

template<typename F>
class FitSineContext : public MinimizeContext
{		
	F& f;	
	const Options& options;
	double _delta_chi_squared = 1.0;
	
public:
	FitSineContext(F& f,const Options& options)
		:f{f},options{options} {	
	}

	double operator() (const std::vector<double>& params) override
	{
		return f({params[0],params[1],params[2]});		
	}
	

	double delta_chi_squared() const override {
		return _delta_chi_squared;
	}
		
	void delta_chi_squared(double delta_chi_squared)
	{
		_delta_chi_squared = delta_chi_squared;
	}
};




template<typename  F>
MinimizerResult fit_sine_Minuit2(std::array<double,3> x0,F chi2f,const Options& options)
{		
	FitSineContext<F> context(chi2f,options);
		
	std::vector<double> p0 {x0[0],x0[1],x0[2]};
	std::ostringstream oss;
	auto result = minimize_locally_Minuit2(p0,context,oss);
	result.x[1] = period_2pi(result.x[1]);
	return result;
}

#endif


template<typename F>
MinimizerResult fit_sine_grad(std::array<double,3> x0,F chi2f,const Options& )
{	
	const int max_iter=100;		
	
	// minimize
	double h = 1e-5;
	double precision = 1e-8;
	auto lm = minimize_locally_a(x0,chi2f,h,max_iter,precision);	
	if (lm.i == max_iter) {
		h *= 0.1;
		auto lmb = minimize_locally_b(x0,chi2f,h,max_iter,precision);
		if (lmb.y <= lm.y)
			lm = lmb;
	}
	
	// set result and evaluate uncertainties
	lm.x[1] = period_2pi(lm.x[1]);
	
	MinimizerResult result;	
	result.valid = lm.i < max_iter;
	result.y = lm.y;
	result.x.push_back(lm.x[0]);
	result.x.push_back(lm.x[1]);
		
	std::array<double,3> temp_params;
	auto ua_pos = fit_parameter_uncertainty(lm.x,lm.y,1.,-100,100,chi2f,1e-6,temp_params,temp_params.at(0));
	auto ua_neg = fit_parameter_uncertainty(lm.x,lm.y,1.,-100,100,chi2f,-1e-6,temp_params,temp_params.at(0));
	auto up_pos = fit_parameter_uncertainty(lm.x,lm.y,1.,-AETHER_PI,AETHER_PI,chi2f,1e-6,temp_params,temp_params.at(1));
	auto up_neg = fit_parameter_uncertainty(lm.x,lm.y,1.,-AETHER_PI,AETHER_PI,chi2f,-1e-6,temp_params,temp_params.at(1));	
	
	auto ua = std::max(ua_pos,ua_neg);
	auto up = std::max(up_pos,up_neg);
	
	result.u.push_back(ua);
	result.u.push_back(up);	
	return result;
	
}


/**
 * \~german
 * Sinus an Daten mit Unsicherheiten anpassen.
 * 
 * \~english
 * Fit sine to data with uncertainties.
 * 
 * \~
 * @param data
 * @param uncertainties
 * @param single
 * @return phase and amplitude
 */
MinimizerResult fit_sine(const std::array<double,17>& data, 
						 const std::array<double,17>& uncertainties,const Options& options)
{
	
	const size_t N = azimuths(options);
		
	// chi squared function object
	auto chi2f = [&data,&uncertainties,&N](const std::array<double,3>& x) {
		double chi=0;
		for (size_t i=0;i<N;i++) {
			// if the sign of the phase x[1] is changed, it changes the residuals statistics
			auto y = x[0]*std::sin(i*AETHER_2PI/8 + x[1])+x[2];
			chi += sqr((data[i]-y)/uncertainties[i]);
		}
		
		return chi;
	};
	
	// find good start point	
	std::array<double,3> x0;	
	auto mm = std::minmax_element(data.begin(),data.end());
	x0[0] = (*mm.second - *mm.first)*0.5; // amplitude
	x0[1] = 0; // phase	
	x0[2] = mean_value(data); // ordinate offset
	
	double min_chi = chi2f(x0);
	double min_phase = x0[1];
	for (int i=1;i<8;i++) {
		x0[1] = i*AETHER_2PI/8;
		auto chi = chi2f(x0);
		if (chi < min_chi) {
			min_chi = chi;
			min_phase = x0[1];
		}			
	}
	x0[1] = min_phase;
	
	
	switch (options.minimizer) {
#ifdef AETHER_MINUIT
	case Options::Minimizer::Minuit2:
		return fit_sine_Minuit2(x0,chi2f,options);		
#endif
	case Options::Minimizer::Grad:
		return fit_sine_grad(x0,chi2f,options);		
	default:				
		throw std::runtime_error("unknown minimizer");
	}		
}



void signal_stats(const MinimizerResult& result,const std::vector<ExtractedSignal>& extracted_signals,
				  const Options& options,std::ostream& os)
{		
	auto theory = create_theory(options.theory);
	TheoryParameters params {result.x.at(0),result.x.at(1),result.x.at(2)};
	SignalStats stats;
	chi_squared_sum(*theory,params,extracted_signals,options,true,stats);
		
	int skipped=0;
	os << std::setprecision(3) << std::fixed;
	os << " #    χ²      Signal\n";	
	for (size_t i=0;i<extracted_signals.size();i++) {
		while (options.disabled_signals.find(i+1+skipped)!=options.disabled_signals.end()) {
			os << std::setw(2) << i+1+skipped << "   disabled\n";		
			skipped++;
		}		
		os << std::setw(2) << i+1+skipped;		
		os << std::setw(8) << stats.chis.at(i);				
		os << std::setw(5) << epoch(extracted_signals.at(i).left_data_sheets.at(0)); // epoch is always the same
		os << " ";
		auto& left_ds = extracted_signals.at(i).left_data_sheets;
		output_separated(os,left_ds.begin(),left_ds.end(),",",[](const DataSheet& data_sheet){
			return data_sheet.no;
		});
		os << " - ";
		auto& right_ds = extracted_signals.at(i).right_data_sheets;
		output_separated(os,right_ds.begin(),right_ds.end(),",",[](const DataSheet& data_sheet){
			return data_sheet.no;
		});		
		os << "\n";	  
				
	}
	
	// dont forget the last ones
	int i = extracted_signals.size();
	while (options.disabled_signals.find(i+1+skipped)!=options.disabled_signals.end()) {
		os << std::setw(2) << i+1+skipped << "   disabled\n";		
		skipped++;
	}
	
	// more stats
	os << "\n";	
	os << std::setw(10) << mean_value(stats.chis.begin(),stats.chis.end()) << "  Mean\n";
	
	// Be carefull! chis are sorted after this!
	std::sort(stats.chis.begin(),stats.chis.end());		
	os << std::setw(10) << median_value(stats.chis) << "  Median\n";
	os << "\n";
		
		
	// Residuals		
	auto meanR = mean_value(stats.residuals.begin(),stats.residuals.end());
	auto varR = sample_variance(stats.residuals.begin(),stats.residuals.end(),meanR);
	os << "Residuals distribution: µ = " << meanR << ", σ² = " << varR << ", ";	
	std::sort(stats.residuals.begin(),stats.residuals.end());			
	double A = test_for_normality(stats.residuals);		
	//os << std::defaultfloat << std::setprecision(6); // p-value precision is low
	os << "p-value = " << p_value_A(A) << "\n";
	if (options.residuals) {
		output_separated(std::cout, stats.residuals, "\n"); 
		std::cout << "\n";
	}	
}



void fit(const std::vector<SignalExtractionExpression>& expressions,
		 const std::vector<DataSheet>& data_sheets,const Options& options)
{	
	std::ostream& os = fit_ostream(options);

	//output_expressions(std::cout,expressions);
	os << "\nExtracting signals..." << std::endl;
	auto extracted_signals = extract_signals(expressions,data_sheets,options);	
	os << "Fitting theory to " << extracted_signals.size() << " signals.\n\n";

	TheoryParameters start_params;
	if (options.start_params) {
		os << "Using starting point:";
		start_params = *options.start_params;
	}
	else {
		os << "Scanning for good starting point..." << std::flush;
		start_params = best_of_random_samples(extracted_signals,options,ProgressInPercent(os));
	}

	os << "\n";
	os << "v = " << start_params.v << " m/s\n";
	os << "α = " << rad_to_h(start_params.a) << " h\n";
	os << "δ = " << deg(start_params.d) << " °\n";
	os << "\n";

	os << "Fitting..." << std::flush;
	MinimizerResult result = fit_theory(start_params,extracted_signals,options);
	os << "\n";	

	write_fit_result(os,result,false,extracted_signals.size(),options);
	os << "\n";
	if (options.stats || options.residuals) {
		signal_stats(result,extracted_signals,options,os);
	}
	
	if (!result.valid) {
		os << "WARNING: Minimizing did not converge\n";
	}

	if (options.contour) {		
		// confidence
		// p     | ∆χ² (1-3 params)
		//--------------------------------
		// 68.3% | 1     2.3    3.5    (1-Sigma)
		// 90%   | 2.7   4.6    6.3
		// 95.5% | 4     6.2    8.0    (2-Sigma)
		// 99%   | 6.6   9.2   11.3
		// 99.7% | 9    11.8   14.2    (3-Sigma)
		if (result.valid) {
			write_fit_result(std::cout,result,true,extracted_signals.size(),options);
			write_coordinates(std::cout,result);
			contour(result,extracted_signals,options,std::cout);
		}
		else {
			std::cerr << "Skipped contour scanning\n";
		}
	}
}



void execute_aggregate_fit(const std::vector<DataSheet>& data_sheets,const Options& options)
{
	if (!options.data_filename.empty() || options.subtract_data) {
		std::cerr << "Data file processing not allowed\n";
		throw ExitException();
	}

	// Normal output is redirected to stderr, if contour data is requested.
	std::ostream& os = fit_ostream(options);

	os << "Enter expressions to extract the signals from the data sheets.\n";
	os << "Type 'fit' to commit and start the fitting process.\n";
	os << "Type 'cancel' to cancel the input.\n";

	int ln = 0;
	std::vector<SignalExtractionExpression> expressions;
	std::string line;
	do {
		os << ": ";
		std::getline(std::cin, line);
		if (!std::cin)
			break;
		ln++;
		if (line=="fit" || line=="cancel")
			break;
		auto pure_line = trim(remove_comments(line));
		if (pure_line.empty())
			continue;

		try {
			auto expr = parse_fit_expression(pure_line);
			if (!validate_fit_expression(os,expr)) {
				os << "Line " << ln << ": " << line << "\n"; 
			}
			expressions.push_back(std::move(expr));
		}
		catch(std::invalid_argument& e) {
			os << "\nExpression ignored: " << e.what() << "\n";
			os << "Line " << ln << ": " << line << "\n";
		}
	} while(true);
	
	if (line == "fit" || std::cin.eof()) {
		if (!expressions.empty()) {
			fit(expressions,data_sheets,options);
		}
	}

}



template<typename F>
std::map<double,double>
aggregate_mean_spectrum(const std::vector<DataSheet>& data_sheets,const Options& options,F dft)
{	
	std::map<double,double> aggregated_amplitude;
	std::map<double,size_t> freq_counts;
	std::map<double,double> current_max_amps;
	
	const int refn = 20; // 20 turns as reference, for frequency bins	
	
	for (const auto& data_sheet : data_sheets) {
		// collect max amplitude of all freq, into ref freq bins
		current_max_amps.clear();
		std::map<double,std::complex<double>> result = dft(data_sheet,options);
		for (const auto& entry : result) {
			auto& f = entry.first;
			auto& z = entry.second;	
			auto reff = std::round(f*refn)/refn; // reference freq bin			
			current_max_amps[reff] = std::max(current_max_amps[reff], std::abs(z)/20); // 1/10 fringes to wave length			
		}		
		for (auto& entry : current_max_amps) {
			auto& f = entry.first;
			auto& A = entry.second;	
			aggregated_amplitude[f] += A;
			freq_counts[f] += 1;
		}
	}
	
	for (const auto& entry : freq_counts) {
		auto& f = entry.first;
		auto& n = entry.second;
		aggregated_amplitude[f] /= n;
	}
	
	aggregated_amplitude.erase(0);
	return aggregated_amplitude;
}



void execute_aggregate_mean_raw_spectrum(const std::vector<DataSheet>& data_sheets,const Options& options)
{	
	auto aggregated_amplitude = aggregate_mean_spectrum(data_sheets,options,raw_spectrum);
		
	std::cout << "# freq\tamplitude\n";
	for (const auto& entry : aggregated_amplitude) {
		std::cout << entry.first << "\t" << entry.second << "\n";
	}	
}



void execute_aggregate_mean_spectrum(const std::vector<DataSheet>& data_sheets,const Options& options)
{		
	auto aggregated_amplitude = aggregate_mean_spectrum(data_sheets,options,spectrum);
			
	std::cout << "# freq\tamplitude\n";
	for (const auto& entry : aggregated_amplitude) {
		std::cout << entry.first << "\t" << entry.second << "\n";
	}	
}



void execute_aggregate_mean(Action action,const std::vector<DataSheet>& data_sheets,const Options& options)
{
	switch(action) {		
	case Action::Reduce:
		execute_aggregate_mean_reduce(data_sheets,options);
		break;
	case Action::Spectrum:
		execute_aggregate_mean_spectrum(data_sheets,options);
		break;
	case Action::RawSpectrum:
		execute_aggregate_mean_raw_spectrum(data_sheets,options);
		break;
	default:
		std::cerr << "Given action can not be aggregated\n";
		throw ExitException();
	}
}







std::mt19937 simulation_rengine;	

/**
 * \~german
 * Simulierte Daten setzen
 * 
 * Die Messdaten eines Datenblattes werden durch simulierte Daten ersetzt.
 * Die simulierten Daten bestehen aus dem theoretischen Signal, einem konstanten
 * systematischen Fehler, und normal verteilten statistischen Fehlern.
 * 
 * Die Fehler sollen:
 * - ähnliches Bild erzeugen
 * - ähnliche Unsicherheit erzeugen
 * - keine systematischen Restfehler im Differenzsignal erzeugen 
 * - ähnliche Ergebnisse beim Test auf Normalverteilung liefern
 * Referenz sind die guten Datenblätter.
 * 
 * Alle Datenblattattribute werden als wahr angenommen 
 * und das Vorzeichen wird entsprechend gesetzt.
 * 
 * \~english
 * 
 * The measured data is replaced with simulated data.
 * The simulated data consists of the theoretical signal, a constant
 * systematic error, and normal distributed statistical errors.
 * 
 * The error should:
 * - create a similar image
 * - create a similar uncertainty
 * - create no systematic error residue in a difference signal
 * - create similar results at a test for normality
 * The good data sheets are used as references.
 * 
 * All data sheet attributes are considered to be true
 * and the sign will be set accordingly.
 * 
 * \~
 * @param data_sheet
 * @param options
 */
void set_simulated_data(DataSheet& data_sheet,const Options& options)
{
	auto theory = create_theory(options.theory);
	auto theory_displs = fringe_displacements(*theory,options.theory_params,data_sheet,options);
	auto z = DFT_analyze(2,theory_displs.begin(),theory_displs.end()-1);
	auto A = std::abs(z);
	auto phi = std::arg(z);
	
	const double sigma_p = 0.8;
	const double sigma_a = 0.5;
	std::normal_distribution<double> distp(0,sigma_p);
	std::normal_distribution<double> dista(0,sigma_a);
		
	const double sigma = 0.015;
	std::normal_distribution<double> dist(0,sigma);	

	const double p_step = options.sim_simple ? 0. : 0.04;
	double p_shift = -double(data_sheet.turns.size())/2. * p_step;
	
	std::array<double,17> displacements = theory_displs;
	for (auto& turn : data_sheet.turns) {				
					
		if (options.sim_sys) {											
			auto p = -AETHER_PI/4; // phase shift -pi/4 = sys error signal source at azimuth 12
			auto a = 0.02;
			auto err_p = options.sim_simple ? 0. : distp(simulation_rengine);
			auto err_a = options.sim_simple ? 0. : dista(simulation_rengine);		
			p += err_p;		
			a = a*(1 + err_a);
			
			set_sine(a,p+p_shift,0,8,displacements); // double period									
			add_sine(a*1.5,p+p_shift+err_p,0,16,displacements); // single period, also creates non normal distribution 
			//add_array(displacements,theory_displs);
			add_sine(A,phi+AETHER_PI_2+p_shift,0,8,displacements);
		}					
		
		auto offs = turn.distances[0];
		auto drift = turn.distances[16]-turn.distances[0]; // needed for rounding
		for (size_t i=0;i<17;i++) {						
			// rounding also creates some non normal distribution
			turn.distances[i] = std::round(offs + drift*i/16. + (displacements[i] + dist(simulation_rengine))*20); 
		}		
		
		int sign = data_sheet.reverse ? -1 : 1;
		if (turn.reverse || turn.invert)
			sign = -sign;
		for (auto& d : turn.distances)
			d *= sign;	
		
		p_shift += p_step;
	}	
	
}



//-----------------------------------------------------------------
//    aggregate signals
//-----------------------------------------------------------------





constexpr bool is_zero_dT(double dT)
{	
	return dT > -0.05 && dT < 0.05;
}



struct SheetDiffStats
{
	double TD = NAN;
	double dT = NAN;	
	bool day_and_night = false;
	bool sun_low = false;
};



template<typename Iter>
void temperature_differences(Iter lbegin,Iter lit,Iter lend,
							 Iter rbegin,Iter rit,Iter rend,Iter pit,SheetDiffStats& diff_stats)
{
	const DataSheet& left = lit->get();
	const DataSheet& right = pit->get();
	
	// TD (simple)
	/*auto left_TD = max_TD(left);	
	auto right_TD = max_TD(right);
	if (left_TD && right_TD) {
		diff_stats.TD = std::abs(*left_TD-*right_TD); 
	}*/
	
	// TD (cross difference)
	auto left_th = mean_thermometers(left.thermometers_start,left.thermometers_end);
	auto right_th = mean_thermometers(right.thermometers_start,right.thermometers_end);
	if (left_th && right_th) {		
		auto ltd = left_th->N + left_th->S - (left_th->W + left_th->E);
		auto rtd = right_th->N + right_th->S - (right_th->W + right_th->E);
		diff_stats.TD = std::abs(ltd-rtd);
	}
		
					
	// dT
	auto left_dT = mean_dT(left);
	if (!left_dT) {
		if (lit==lbegin && lit+1==lend)
			return;
		int offs = lit==lbegin ? 1 : -1; // check previous			
		left_dT = mean_dT(*(lit+offs),*lit);
		if (std::isnan(*left_dT))
			return;
	}
	
	auto right_dT = mean_dT(right);
	if (!right_dT) {
		if (rit==rbegin && rit+1==rend)
			return;
		int offs = rit==rbegin ? 1 : -1; // check previous		
		right_dT = mean_dT(*(rit+offs),*rit);
		if (std::isnan(*right_dT))
			return;
	}
	
	if (!(is_zero_dT(*left_dT) || is_zero_dT(*right_dT))) {
		if ((*left_dT<0 && *right_dT>0) || (*left_dT>0 && *right_dT<0))
			return; // T change in different directions is not allowed		
	}
	
	diff_stats.dT = std::abs(*left_dT-*right_dT);
			
}



void time_of_day_diff(const DataSheet& prev, const DataSheet& next,SheetDiffStats& diff_stats)
{
	auto tod1 = time_of_day(prev);
	auto tod2 = time_of_day(next);			
	diff_stats.sun_low = 
			tod1==TimeOfDay::Sunset || tod1==TimeOfDay::Sunrise || 
			tod2==TimeOfDay::Sunset || tod2==TimeOfDay::Sunrise;
	
	diff_stats.day_and_night = (tod1==TimeOfDay::Night) != (tod2==TimeOfDay::Night);	
}



bool selected(const SheetDiffStats& diff_stats,const Options& options)
{	
	// uses round() to match displayed values in statistics etc.
	
	if (std::isnan(diff_stats.TD))
		return false;
	if (std::round(diff_stats.TD*100) > std::round(options.signals_dTD*100)) // 1/100 °C
		return false;
	
	if (std::isnan(diff_stats.dT))
		return false;
	if (std::round(diff_stats.dT*100) > std::round(options.signals_ddT*100)) // 1/100 °C
		return false;	
	
	if (!options.low_sun && diff_stats.sun_low)
		return false;
	if (!options.day_and_night && diff_stats.day_and_night)
		return false;
			
		
	return true;
}



struct AggregatedDiffStats
{	
	double TD = INFINITY; // init for min()
	double dT = INFINITY;	
	int n = 0;
};



/**
 *
 * \~
 * \warning O(n)=n!
 */
template<typename Iter>
bool selected(Iter lbegin,Iter lfirst,Iter llast,Iter lend,
			  Iter rbegin,Iter rfirst,Iter rlast,Iter rend,const Options& options,
			  AggregatedDiffStats& aggregated_diff_stats)
{
	std::vector<std::reference_wrapper<const DataSheet>> pbuffer(rfirst,rlast+1); // is sorted in call context
	std::vector<SheetDiffStats> diffs;
	
	// Compare the sequence of left sheets with all permutations of the right sheets.
	// If all sheets of one permutation are selected, that is a success.
	do {
		diffs.clear();		
		Iter lit = lfirst;		
		Iter pit = pbuffer.begin();
		while (lit!=llast+1 && pit!=pbuffer.end()) {			
			auto rit = std::find_if(rfirst,rlast+1,[&pit](const DataSheet& data_sheet) {
				return data_sheet.no == pit->get().no;
			});
			SheetDiffStats diff_stats;
			temperature_differences(lbegin,lit,lend,rbegin,rit,rend,pit,diff_stats);
			time_of_day_diff(*lit,*rit,diff_stats);
			if (!selected(diff_stats,options))				 
				break;			
			
			diffs.push_back(std::move(diff_stats));
			++lit;			
			++pit;
		}
			
		if (diffs.size() == size_t(llast-lfirst+1)) { // all selected?			
			auto max_TD_stats = std::max_element(diffs.begin(),diffs.end(),[](const SheetDiffStats& a,const SheetDiffStats& b){
				return a.TD<b.TD;
			});
			auto max_dT_stats = std::max_element(diffs.begin(),diffs.end(),[](const SheetDiffStats& a,const SheetDiffStats& b){
				return a.dT<b.dT;
			});
			aggregated_diff_stats.TD = max_TD_stats->TD;
			aggregated_diff_stats.dT = max_dT_stats->dT;
			aggregated_diff_stats.n = diffs.size();
		
			return true;
		}
	} 
	while (std::next_permutation(pbuffer.begin(),pbuffer.end()));	
	
	return false;
}



bool better_than(const AggregatedDiffStats& a, const AggregatedDiffStats& b)
{
	if (a.n>b.n) {	
		return true;
	}
	
	if (a.n==b.n) {
		// yes, these are different physical quantities
		if (a.TD+a.dT < b.TD+b.dT) {
			return true;
		}
	}
	
	return false;
}



/**
 * \~german
 * Sternzeitlicher Abstand zweier Datenblätter
 * 
 * \~english
 * 
 * \~
 * @param a
 * @param b
 * @return distance in h
 */
double sidereal_distance(const DataSheet& a,const DataSheet& b)
{	
	// longitude is irrelevant in this case
	auto theta1 = sidereal_time(calendar_date(a),0);
	auto theta2 = sidereal_time(calendar_date(b),0);
	return rad_to_h(periodic_distance(theta1,theta2,AETHER_2PI));
}



std::vector<SignalExtractionExpression>
group_signals(const std::vector<std::reference_wrapper<const DataSheet>>& group,const Options& options)
{	
	std::vector<SignalExtractionExpression> signals;	
	if (group.size()<2)
		return signals;
			
	AggregatedDiffStats best_diff_stats;
		
	for (auto it1 = group.begin();it1!=group.end()-1;++it1) {
		for (auto it2 = it1+1;it2!=group.end();++it2) {			
			if (sidereal_distance(*it1,*it2) < options.signals_dt) // minimum distance in time
				continue; // No max distance check needed, largest groups span about 6 h.
			
			int offs = 0; // used to build sub sequences
			do {
				AggregatedDiffStats diff_stats;
				if (selected(group.begin(),it1-offs,it1,group.end(), group.begin(),it2,it2+offs,group.end(),options,diff_stats)) {
					if (better_than(diff_stats,best_diff_stats) && offs>0) { // at least two
						best_diff_stats = diff_stats;
								
						SignalExtractionExpression signal;
						signal.epoch = month_interval(*it1); // same for all
						signal.left.push_back({(it1-offs)->get().no,(it1)->get().no});
						signal.right.push_back({(it2)->get().no,(it2+offs)->get().no});
						if (signals.empty()) 
							signals.push_back(signal);
						else {
							signals.back() = signal; // just one signal currently
						}						
					}
				}
												
				offs++;				
			}
			while(offs <= it1-group.begin() && offs < group.end()-it2);
		}
	}
	
	
	return signals;
}



bool selected(const std::vector<std::reference_wrapper<const DataSheet>>& group1,
			  const std::vector<std::reference_wrapper<const DataSheet>>& group2)
{
	auto sheet1 = group1.front().get();
	auto sheet2 = group2.front().get();
	return epoch(sheet1)==epoch(sheet2) && sheet1.desk_in_sw==sheet2.desk_in_sw;
}



std::vector<SignalExtractionExpression>
intergroup_signals(const std::vector<std::reference_wrapper<const DataSheet>>& group1,
				   const std::vector<std::reference_wrapper<const DataSheet>>& group2,
				   const Options& options)
{
	std::vector<SignalExtractionExpression> signals;	
	AggregatedDiffStats best_diff_stats;
		
	for (auto it1 = group1.begin();it1!=group1.end();++it1) {
		for (auto it2 = group2.begin();it2!=group2.end();++it2) {						
						
			int offs = 0; // used to build sub sequences in groups			
			do {
				auto dtheta = sidereal_distance(*(it1-offs),*(it2+offs)); 
				if (dtheta < options.signals_dt || dtheta > 8.0) 
					break;
				
				AggregatedDiffStats diff_stats;
				if (selected(group1.begin(),it1-offs,it1,group1.end(), group2.begin(),it2,it2+offs,group2.end(),options,diff_stats)) {
					if (better_than(diff_stats,best_diff_stats) && offs>0) { // at least two
						best_diff_stats = diff_stats;
								
						SignalExtractionExpression signal;
						signal.epoch = month_interval(*it1); // same for all
						signal.left.push_back({(it1-offs)->get().no,(it1)->get().no});
						signal.right.push_back({(it2)->get().no,(it2+offs)->get().no});
						if (signals.empty()) 
							signals.push_back(signal);
						else {
							signals.back() = signal; // just one signal currently
						}						
					}
				}
								
				offs++; 
			}
			while(offs <= it1-group1.begin() && offs < group2.end()-it2);
		}
	}
		
	return signals;
}





std::vector<std::vector<SignalExtractionExpression>>
intergroup_signals(const std::vector<std::vector<std::reference_wrapper<const DataSheet>>>& groups,const Options& options)
{
	std::vector<std::vector<SignalExtractionExpression>> signals;
	
	if (groups.size()<2)
		return signals;
	
	// for all combinations of two groups
	for (auto it1 = groups.begin();it1!=groups.end()-1;++it1) {
		for (auto it2 = it1+1;it2!=groups.end();++it2) {
			if (selected(*it1,*it2)) {				
				auto signals12 = intergroup_signals(*it1,*it2,options);
				if (!signals12.empty()) {
					signals.push_back(signals12);
				}
			}
		}
	}
		
	return signals;
}





bool same_group(const DataSheet& previous, const DataSheet& current)
{
	const double max_dt = 1.5; // in h
	
	auto jd1 = julian_date(calendar_date(previous));	
	auto jd2 = julian_date(calendar_date(current));
	if ((jd2-jd1)*24 > max_dt)
		return false;
	return epoch(previous)==epoch(current) && previous.desk_in_sw==current.desk_in_sw;
}



std::vector<std::vector<std::reference_wrapper<const DataSheet>>> 
split_into_groups(const std::vector<DataSheet>& data_sheets)
{
	std::vector<std::vector<std::reference_wrapper<const DataSheet>>> groups;
	
	groups.emplace_back();
	for (const auto& data_sheet : data_sheets) {
		if (!groups.back().empty() && !same_group(groups.back().back(),data_sheet)) {
			groups.emplace_back();
		}
		groups.back().push_back(std::ref(data_sheet));
	}
	
	return groups;
}



void execute_aggregate_signals(std::vector<DataSheet>& data_sheets,const Options& options)
{
	std::sort(data_sheets.begin(),data_sheets.end());
	
	auto groups = split_into_groups(data_sheets);
			
	/*std::cerr << "found groups: ";
	for (const auto& group : groups) {
		 std::cerr << epoch(group.front().get()) << "-" << group.front().get().no << " ";
	}
	std::cerr << "\n";*/
	
	std::cout << "##################################################\n";
	std::cout << "# Generated signal extractions using options\n";
	std::cout << "# -signals_dTD   " << options.signals_dTD << " °C\n";
	std::cout << "# -signals_ddT   " << options.signals_ddT << " °C per ¼h\n";
	std::cout << "# -signals_dt    " << options.signals_dt << " h\n";
	std::cout << "# -day_and_night " << yesno(options.day_and_night) << "\n";
	std::cout << "# -low_sun       " << yesno(options.low_sun) << "\n";
	std::cout << "##################################################\n";
	
	std::cout << "\n";
		
	std::cout << "### Group signals ###\n\n";
	for (const auto& group : groups) {		
		auto signals = group_signals(group,options);		
		if (!signals.empty()) {	
			output_expressions(std::cout,signals);						
		}
	}
	
	std::cout << "\n";
	
	std::cout << "### Intergroup signals ###\n\n";
	auto groups_signals = intergroup_signals(groups,options);
	for (const auto& signals : groups_signals) {
		if (!signals.empty()) {			
			output_expressions(std::cout,signals);						
		}
	}
	
	std::cout << "\n";
}



void execute(Action action,const Options& options,std::vector<DataSheet>& data_sheets)
{
	switch(options.aggregation_method) {
	case Options::AggregationMethod::List:
		execute_aggregate_list(data_sheets,options);
		break;
	case Options::AggregationMethod::Test:
		execute_aggregate_test(data_sheets,options);
		break;
	case Options::AggregationMethod::Mean:
		execute_aggregate_mean(action,data_sheets,options);
		break;
	case Options::AggregationMethod::Sidereal:
		execute_aggregate_sidereal(data_sheets,options);
		break;
	case Options::AggregationMethod::Diff:
		execute_aggregate_diff(data_sheets,options);
		break;
	case Options::AggregationMethod::Params:
		execute_aggregate_params(data_sheets,options);
		break;
	case Options::AggregationMethod::ModelChi:
		execute_aggregate_model_chi(data_sheets,options);
		break;
	case Options::AggregationMethod::Fit:
		execute_aggregate_fit(data_sheets,options);
		break;		
	case Options::AggregationMethod::Signals:
		execute_aggregate_signals(data_sheets,options);
		break;
	default:
		throw std::runtime_error("unknown aggregation method");
	}

}

}//aether
