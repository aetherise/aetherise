/**
 * \file
 *
 * \~german
 * @brief Datenformate erzeugen oder Daten formatiert ausgeben
 *
 * \~english
 * @brief Generate data formats or output data formatted
 *
 */

#ifndef GENERATORS_H
#define GENERATORS_H

#include "aetherise.h"
#include "DataSheet.h"
#include "data_reduction.h"
#include "models.h"

#include <array>
#include <map>
#include <ostream>

namespace aether {




/**
 * \~german
 * Allgemeine Statistiken eines Datenblattes
 * die nicht von einer bestimmten Aktion abhängen.
 *
 * \~english
 * Common stats of a data sheet
 * that do not depend on an specific action.
 *
 * \~
 * @param os
 * @param data_sheet
 */
void write_data_sheet_stats(std::ostream& os,const DataSheet& data_sheet);

void write_header(std::ostream& os,const DataSheet& data_sheet);

void write_raw(std::ostream& os,const DataSheet& data_sheet);

void write_reduce_raw(std::ostream& os,const std::vector<ReducedTurn>& turns);

void write_reduced_data(std::ostream& os,const DisplacementData& distance_data,
						const DataSheet& data_sheet,const Options& options);


void write_aggregated_data(std::ostream& os,  const DisplacementData& aggregated_distances,
						   const Options& options);


void write_aggregated_data(std::ostream& os,const std::map<int,SiderealData>& aggregated_sidereal_data,
						   const Options& options);


void write_aggregated_data(std::ostream& os,const std::vector<ErrorModelParameters>& params);



/**
 *
 * @param os
 * @param data_sheets
 * @param options
 */
void write_list(std::ostream& os,const std::vector<DataSheet>& data_sheets,const Options& options);



/**
 * \~german
 * Chi-Quadrat-Statistik in gegebenen Ausgabestrom schreiben.
 *
 * \~english
 * Write chi squared statistics to given output stream.
 *
 * \~
 * @param os output stream
 * @param chi χ² value
 * @param f degrees of freedom
 * @param p p-value
 * @param prefix prefix to output in front of each line
 */
void write_chi_squared_stats(std::ostream& os,double chi,int f,double p,const std::string& prefix = {});

}//


#endif // GENERATORS_H
