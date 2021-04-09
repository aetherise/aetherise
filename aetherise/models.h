/**
 * \file
 *
 * \~german
 * @brief Fehlermodelle
 *
 * \~english
 * @brief Error models
 *
 */


#ifndef MODELS_H
#define MODELS_H

#include "DataSheet.h"
#include "utils.h"

#include <array>


namespace aether {


struct ErrorModelParameters
{
	double desk_azimuth; ///< [1,16]
	double desk;	
	double temperature;
};



/**
 * 
 * \~
 * \warning experimental, incomplete, not up to date 
 * 
 * @param data_sheet
 * @param options
 * @return
 */
std::array<double,17>
systematic_error_displacements(const DataSheet& data_sheet,const Options& options);



/**
 *
 * \~
 * \warning experimental, incomplete, not up to date 
 * 
 * @param start
 * @param end
 * @param params
 * @param options
 * @param desk_disabled
 * @return
 */
std::array<double,17>
systematic_error_displacements(const optional<DataSheet::Thermometers>& start,
							  const optional<DataSheet::Thermometers>& end,
							  const ErrorModelParameters& params,
							  const Options& options,bool desk_disabled=false);

}//

#endif // MODELS_H
