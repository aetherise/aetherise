/**
 * \file
 *
 * \~german
 * @brief Die verschiedenen Datenanalyseverfahren
 *
 * \~english
 * @brief The different methods of data analysis
 *
 */

#ifndef AETHERISE_H
#define AETHERISE_H

#include "DataSheet.h"
#include "utils.h"

#include <vector>
#include <string>
#include <random>


/**
 * \~german
 * Der Namensraum für das aetherise Projekt.
 *
 * \~english
 * The namespace for the aetherise project.
 */
namespace aether {

class MinimizerResult;
class Theory;


const double Sidereal_Aggregation_Bin_Width = 0.5; ///< in hours

extern std::mt19937 simulation_rengine;	


/**
 * \~german
 * Streifenverschiebungen der Daten (mit Unsicherheiten) und
 * der Theorie und des Modells.
 *
 * \~english
 * Fringe displacements of the data (with uncertainties) and
 * the theory and the model.
 */
struct DisplacementData {
	std::array<double,17> data;
	std::array<double,17> uncertainties;
	std::array<double,17> theory;
	std::array<double,17> model;
};



/**
 * \~german
 * Eine einzelne reduzierte Messung (Umdrehung)
 *
 * \~english
 * A single reduced measurement (turn)
 *
 */
struct ReducedTurn
{
	int i;
	DataSheet::Turn turn;
	std::array<double,17> displacements;
};



/**
 * \~german
 * Datensatz der bei der Aggregierung nach Sternzeit anfällt.
 *
 * \~english
 * Record that is used for the aggregation by sidereal time.
 */
struct SiderealData
{
	double data_amplitude = 0;
	double theory_amplitude = 0;
	double model_amplitude = 0;
	int n = 1;	
	double TD = 0;
	int nT = 1;
};





void subtract_data(DisplacementData& displacements,const Options& options);
MinimizerResult fit_sine(const std::array<double,17>& data, 
						 const std::array<double,17>& uncertainties,const Options& options);
TheoryParameters add_earth_orbit(const TheoryParameters& params,const DataSheet& data_sheet);
double humidity(const DataSheet& data_sheet);
double index_of_refraction(const DataSheet& data_sheet);
std::array<double,17> fringe_displacements(const Theory& theory,const TheoryParameters& params,
										   const DataSheet& data_sheet,const Options& options);
ConfidenceInterval rejection_quota(int r,int n);
double estimate_amplitude(const std::array<double,17>& data);
int sidereal_time_to_key(double sidereal_time);
double key_to_mean_sidereal_time(int key);


int degrees_of_freedom(int n,const Options& options);

void set_simulated_data(DataSheet& data_sheet,const Options& options);

void execute(Action action,const Options& options,const DataSheet& data_sheet,const std::string& filename);
void execute(Action action,const Options& options,std::vector<DataSheet>& data_sheets);





}//aether

#endif // AETHERDL_H
