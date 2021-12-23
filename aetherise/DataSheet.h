/**
 * \file
 *
 * \~german
 * @brief Datenblatt
 *
 * \~english
 * @brief data sheet
 *
 */

#ifndef DATASHEET_H
#define DATASHEET_H

#include "utils.h"

#include <vector>
#include <array>
#include <string>


namespace aether {




/**
 * \~german
 * Ein Datenblatt von Miller
 *
 * \~english
 * A data sheet from Miller
 */
struct DataSheet
{	
	struct Date {
		int year;
		short month;
		short day;
	};

	struct Time {
		short hour; // 24h
		short minute;
	};

	struct Thermometers {		
		double N,E,S,W; // (°C)
	};

	enum class WeatherObscurations {
		Clear,Haze,Mist,Fog
	};

	struct Turn {
		std::array<float,17> distances; // in 1/10 of a fringe from reference point
		bool adjust = false;
		bool reverse = false;
		bool bad = false;
		bool cancel = false;
		bool invert = false;
		bool plus = false;
		bool minus = false;
	};

	Date date;
	int no;

	optional<Time> thermometers_start_time;
	optional<Time> thermometers_end_time;
	optional<Thermometers> thermometers_start;
	optional<Thermometers> thermometers_end;

	Time mean_observation_time;
	optional<Time> local_mean_time;
	Time sidereal_time;
	optional<double> weight;
	optional<int> fringes;
	bool sign_correct;
	std::string comment;	

	Time start_time;
	Time end_time;
	bool inverted = false;
	bool reverse = false;
	bool desk_in_sw = false;
	bool desk_disabled = false;
	bool visitors = false;
	bool corrugated_paper = false;

	bool weather_wind = false;
	bool weather_rain = false;
	WeatherObscurations weather_obscuration = WeatherObscurations::Clear;
	bool weather_clouds = false;

	std::vector<Turn> turns;
};



/**
 * \~german
 * Statistiken eines Datenblattes
 *
 * \~english
 * Stats of a data sheet
 */
struct DataSheetStats
{
	int number_of_adjust = 0;
	int number_of_plus = 0;
	int number_of_minus = 0;
	double drift = 0.0;
	double abs_drift = 0.0;
	optional<double> mean_T;
	optional<double> mean_dT;
	optional<double> max_TD;
};


/**
 * \~german
 * Die vier Epochen in denen Messungen durchgeführt wurden.
 * 
 * \~english
 * The four epochs in which measurements took place.
 *  
 */
enum class Epoch
{
	Apr,Aug,Sep,Feb
};



/**
 * \~german
 * Vergleicht Datum und Nummerierung
 *
 * \~english
 * Compares date and numbering
 *
 * \~
 * @param a
 * @param b
 * @return a earlier than b?
 */
bool operator <(const DataSheet& a, const DataSheet& b);


bool operator == (const DataSheet::Date& a, const DataSheet::Date& b);
bool operator == (const DataSheet::Time& a, const DataSheet::Time& b);

bool operator <(const DataSheet::Date& a, const DataSheet::Date& b);
bool operator <(const DataSheet::Time& a, const DataSheet::Time& b);

std::ostream& operator<<(std::ostream& os,const DataSheet::Date& date);
std::ostream& operator<<(std::ostream& os,const DataSheet::Time& time);

DataSheet::Thermometers operator -(const DataSheet::Thermometers& a, const DataSheet::Thermometers& b);


/**
 * \~german
 * Statistiken zu einem Datenblatt erzeugen.
 *
 * \~english
 * Create stats for a data sheet.
 *
 * \~
 * @param data_sheet
 * @return
 */
DataSheetStats data_sheet_stats(const DataSheet& data_sheet);



/**
 * \~german
 * Monat zu Epoche wandeln.
 *
 * \~english
 * Transforms month to epoch.
 *
 * \~
 * @param data_sheet
 * @return
 */
Epoch epoch(const DataSheet& data_sheet);
Epoch epoch(short month);

std::ostream& operator << (std::ostream& os,const Epoch& epoch);



/**
 * \~
 * @param data_sheet
 * @return 
 */
IntegerInterval month_interval(const DataSheet& data_sheet);





/**
 * \~german
 * Ein Datenblatt aus Datei laden.
 *
 *
 * \~english
 * Load a data sheet from a file.
 *
 *
 * \~
 * @param filename
 * @return
 */
DataSheet load_data_sheet_csv(const std::string& filename,const Options& options);


std::string data_sheet_attributes(const DataSheet& data_sheet);
std::string data_sheet_weather(const DataSheet& data_sheet);

/**
 * \~german
 * Datenblatt im CSV-Format ausgeben.
 *
 * \~english
 * Output data sheet in CSV format.
 *
 * \~
 * @param os
 * @param data_sheet
 */
void write_data_sheet_csv(std::ostream& os,const DataSheet& data_sheet);



/**
 * \~german
 * Zeit zu Stundenwert im Intervall [0,24) wandeln.
 *
 * \~english
 * Convert Time to a hour value in the interval [0,24).
 *
 * \~
 * @param t
 * @return
 */
double time_to_h(const DataSheet::Time& t);


/**
 * \~german
 * Stunden zu Uhrzeit (24h) wandeln.
 *
 * \~english
 * Convert hours to Time (24h).
 *
 * \~
 * @param h
 * @return
 */
DataSheet::Time h_to_time(double h);


/**
 * \~german
 * Höchste Temperatur der vier Thermometer finden.
 *
 * \~english
 * Find highest temperature of the four thermometers.
 *
 * \~
 * @param th
 * @return
 */
double max_T(const DataSheet::Thermometers& th);



/**
 * \~german
 * Niedrigste Temperatur der vier Thermometer finden.
 *
 * \~english
 * Find lowest temperature of the four thermometers.
 *
 * \~
 * @param th
 * @return
 */
double min_T(const DataSheet::Thermometers& th);



/**
 * \~german
 * Mittlere Temperatur der vier Thermometer.
 *
 * \~english
 * Mean temperature of the four Thermometers.
 *
 * \~
 * @param th
 * @return
 */
double mean_T(const DataSheet::Thermometers& th);



/**
 * \~german
 * Änderung der mittleren Temperatur.
 *
 * \~english
 * Change of the mean temperature.
 *
 * \~
 * @param th1
 * @param th2
 * @return
 */
double mean_dT(const DataSheet::Thermometers& th1,const DataSheet::Thermometers& th2);



/**
 * \~german Hat dT?
 * 
 * Hat das Datenblatt Temperaturablesungen vom Anfang und Ende?
 * 
 * \~english
 * Are there both temperature readings from start and end?
 * 
 * \~
 * @param data_sheet
 * @return 
 */
bool has_dT(const DataSheet& data_sheet);



/**
 * \~german
 * Mittlere Temperatur aller Thermometer.
 *
 * \~english
 * Mean temperature of all thermometers.
 *
 * \~
 * @param data_sheet
 * @return
 */
optional<double> mean_T(const DataSheet& data_sheet);


/**
 * \~german
 * Größter Temperaturunterschied
 * 
 * \~english
 * Highest temperatur difference
 * 
 * \~
 * @param data_sheet
 * @return 
 */
optional<double> max_TD(const DataSheet& data_sheet);



/**
 * \~german
 * Mittlerer Temperaturunterschied
 * 
 * \~english
 * Mean temperature difference
 * 
 * \~
 * @param th1
 * @param th2
 * @return 
 */
optional<double> mean_TD(const optional<DataSheet::Thermometers>& th1,
						 const optional<DataSheet::Thermometers>& th2);


/**
 * \~german
 * Mittlere Temperaturänderung. Normiert auf 15 min.
 * 
 * \~english
 * Mean change of temperature. Normalized to 15 min.
 * 
 * \~
 * @param data_sheet
 * @return 
 */
optional<double> mean_dT(const DataSheet& data_sheet);



/**
 * \~german
 * Mittlere Temperaturänderung zwischen zwei Datenblättern. 
 * Normiert auf 15 min.
 * 
 * \~english
 * Mean change of temperature between two data sheets. 
 * Normalized to 15 min.
 * 
 * \~
 * @param prev
 * @param next
 * @return 
 */
double mean_dT(const DataSheet& prev,const DataSheet& next);



/**
 *
 * \~
 * @param start
 * @param end
 * @return 
 */
optional<DataSheet::Thermometers>
mean_thermometers(const optional<DataSheet::Thermometers>& start, 
				  const optional<DataSheet::Thermometers>& end);



/**
 * \~german
 * Gültigskeitsprüfung der Datenblätter.
 * Meldet Fehler und warnt vor möglichen Fehlern.
 *
 * \~english
 * Validation of the data sheets.
 * Reports errors and warnings of possible errors.
 *
 * \~
 * @param data_sheet
 * @param filename
 * @param os
 * @return
 */
bool validate(const DataSheet& data_sheet,const std::string& filename,std::ostream& os);



/**
 * \~german
 * Tageszeit auf dem Mount Wilson
 *
 * \~english
 * Time of day on the Mount Wilson
 *
 * \~
 * @param data_sheet
 * @return
 */
TimeOfDay time_of_day(const DataSheet& data_sheet);



/**
 * \~german
 * Kalendardatum
 *
 * Datum und Uhrzeit in Kalender wandeln.
 *
 * \~english
 * Transform date and time to a calendar date.
 *
 * \~
 * @param data_sheet
 * @param lmt use LMT from data sheet
 * @return
 */
Calendar calendar_date(const DataSheet& data_sheet,bool lmt=false);




/**
 * \~german 
 * Ist das Datenblatt aus Cleveland ab 1927?
 * 
 * \~english 
 * Is the data sheet from Cleveland since 1927?
 * 
 * \~
 * @param data_sheet
 * @return 
 */
bool is_from_Cleveland(const DataSheet& data_sheet);

}// aether

#endif // DATASHEET_H
