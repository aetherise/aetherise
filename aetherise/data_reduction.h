/**
 * \file
 *
 * \~german
 * @brief Methoden zur Datenreduzierung
 *
 *
 * \~english
 * @brief Methods for data reduction
 *
 */


#ifndef DATA_REDUCTION_H
#define DATA_REDUCTION_H

#include "DataSheet.h"
#include "utils.h"

#include <array>
#include <memory>

namespace aether {



/**
 * \~german
 * Anzahl der Azimute die nach der Reduzierung übrig bleiben.
 * 
 * \~english
 * Number of azimuths that remain after data reduction.
 * 
 * \~
 * @param options
 * @return 
 */
constexpr int azimuths(const Options& options)
{
	return options.single ? 8 : 16;
}




/**
 * \~german
 * Ausgewählte und transformierte Umdrehungen
 *
 * Die Optionen und Filter auf die einzelnen Umdrehungen des Datenblatts anwenden.
 *
 * \~english
 * Apply options and filter to the turns of a data sheet.
 *
 * \~
 * @tparam F Funktor. Will be called with (rownumber, original turn, modified distances).
 */
template<typename F>
void selected_and_transformed_turns(const DataSheet& data_sheet,const Options& options,F f)
{
	short int sheet_sign = 1;
	if (data_sheet.reverse)
		sheet_sign = -sheet_sign;
	if (options.invert_data)
		sheet_sign = -sheet_sign;

	int i = 0;
	for (const DataSheet::Turn& turn : data_sheet.turns) {
		i++;
		if (turn.cancel || turn.bad)
			continue;

		short int turn_sign = sheet_sign;
		if (turn.reverse || turn.invert)
			turn_sign = -turn_sign;

		auto distances = turn.distances; // work with copy
		for (auto& d : distances)
			d *= turn_sign;

		f(i,turn,distances);
	}
}



/**
 * \~german
 * Reduzierte Daten mit Standardmessunsicherheiten
 *
 * \~english
 * Reduced data with standard uncertainties
 *
 */
struct ReducedData
{	
	std::array<double,17> displacements;
	std::array<double,17> uncertainties;
};



/**
 * \~german
 * Methode der Datenreduzierung.
 *
 * Die Daten eines Datenblattes werden von bekannten systematischen Fehlern befreit
 * und zu einem Datensatz zusammengefasst.
 * Unsicherheiten werden ermittelt.
 *
 * \~english
 * The method to reduce the data of a data sheet to a single record.
 * Uncertainties are evaluated.
 */
class ReductionMethod
{
public:
	virtual ~ReductionMethod() = default;

	virtual ReducedData reduce(const DataSheet& data_sheet, const Options& options) const = 0;
};



/**
 * \~german
 * Millers Methode
 *
 * \~english
 * Miller's method
 */
class MillersReduction : public ReductionMethod
{
public:
	ReducedData reduce(const DataSheet &data_sheet, const Options &options) const override;
};



/**
 * \~german
 * Jede Umdrehung erst normieren und dann zusammenfassen.
 * Diese Methode ist gleichwertig mit Millers.
 *
 * \~english
 * Each turn is normalised first and further reduced second.
 * This method is equivalent to Miller's.
 */
class SeparateReduction : public ReductionMethod
{
public:
	ReducedData reduce(const DataSheet &data_sheet, const Options &options) const override;
};



/**
 * No a real reduction Method. 
 * Attempt to replicate Roberts' model of the "systematic drift".
 */
class Roberts2006 : public ReductionMethod
{
public:
	ReducedData reduce(const DataSheet &data_sheet, const Options &options) const override;
};



/**
 * \~german
 * Instanz der gewählten Methode erzeugen.
 *
 * \~english
 * Create an instance of the choosen method.
 *
 * \~
 * @param method
 * @return
 */
std::unique_ptr<ReductionMethod> create_reduction_method(Options::DataReductionMethod method);



/**
 * \~german
 * Datenreduzierung
 *
 * Reduziert die Daten des Datenblattes mit der über die Optionen gewählten Methode.
 *
 * \~english
 * Reduces the data of the data sheet with the method selected in the options.
 *
 * \~
 * @param data_sheet
 * @param options
 * @return
 */
ReducedData reduce_data(const DataSheet& data_sheet,const Options& options);



/**
 * \~german
 * Drift und Versatz von Messdaten entfernen.
 *
 * \~english
 * Remove drift and offset.
 *
 * \~
 * @param distances
 * @return
 */
std::array<double,17> reduce_from_drift_and_offset(const std::array<float,17>& distances);



/**
 *\~german
 * Doppelperiode weiter auf Einzelperiode reduzieren.
 *
 * \~english
 * Reduce double period to single period.
 *
 * \~
 * @param displacements
 */
void reduce_to_single_period(std::array<double,17>& displacements);



/**
 * \~german
 * Doppelperiode weiter auf Einzelperiode reduzieren,
 * mit Fehlerfortpflanzung.
 *
 * \~english
 * Reduce double period to single period,
 * with uncertainty propagation.
 *
 * \~
 * @param reduced_data
 */
void reduce_to_single_period(ReducedData& reduced_data);


}//aether

#endif // DATA_REDUCTION_H
