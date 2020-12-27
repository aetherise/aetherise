/**
 * \file
 *
 * \~german
 * @brief Kommandozeile auswerten
 *
 * \~english
 * @brief Command line parsing
 *
 */

#ifndef CMD_LINE_H
#define CMD_LINE_H

#include "utils.h"
#include "Filter.h"

#include <string>
#include <vector>

namespace aether {




struct SignalExtractionExpression
{
	IntegerInterval epoch;
	std::vector<IntegerInterval> left,right;
};


/**
 * \~german
 * Argumente der Kommandozeile einlesen.
 *
 * \~english
 * Parses the arguments of the command line.
 *
 * \~
 * @param[in] argc
 * @param[in] argv
 * @param[out] filter
 * @param[out] action
 * @param[out] aggregate
 * @param[out] options
 * @param[out] filenames
 */
void parse_arguments(int argc,char* argv[],Filter& filter,Action& action,bool& aggregate,Options& options,
					 std::vector<std::string>& filenames);

/**
 *
 * @param line
 * @return
 */
SignalExtractionExpression parse_fit_expression(const std::string& line);

}//


#endif // CMD_LINE_H
