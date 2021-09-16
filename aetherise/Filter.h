/**
 * \file
 *
 * \~german
 * @brief Datenblätter filtern
 *
 * \~english
 * @brief Filter data sheets
 *
 */


#ifndef FILTER_H
#define FILTER_H


#include <limits>
#include <vector>


namespace aether {

class DataSheet;
class Options;



/**
 * \~german
 * Ein Filter für Datenblätter
 *
 * \~english
 * A filter for data sheets
 */
struct Filter
{
	/**
	 * \~german Intervall
	 * \~english Interval
	 */
	struct Interval {
		double min = -std::numeric_limits<double>::infinity();
		double max = std::numeric_limits<double>::infinity();

		Interval()=default;
		Interval(double min,double max)
			:min{min},max{max} {
		}
	};

	std::vector<Interval> no;
	std::vector<Interval> year,month,day;
	std::vector<Interval> T; // temperature
	std::vector<Interval> TD; // temperature difference
	std::vector<Interval> dT; // temperature change
	std::vector<Interval> mean_dT;
	std::vector<Interval> adjust;
	std::vector<Interval> time;
	std::vector<Interval> sidereal_time;
	std::vector<Interval> weight;
	std::vector<Interval> fringes;

	std::vector<Interval> amplitude;
	std::vector<Interval> drift;
	std::vector<Interval> abs_drift;
	std::vector<Interval> uncertainty;
	std::vector<Interval> theory_amplitude;
	bool nw = false;
	bool sw = false;
    bool sign_correct = false;
    bool sign_correct_missing = false;
};



/**
 * \~german
 * Ist das Datenblatt von dem gegebenen Filter ausgewählt?
 *
 * \~english
 * Is the data sheet selected by the given filter?
 *
 * \~
 * @param data_sheet
 * @param options
 * @param filter
 * @return
 */
bool selected(const DataSheet& data_sheet,const Options& options,const Filter& filter);



}//aether

#endif // FILTER_H
