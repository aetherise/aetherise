/**
 * \file
 *
 * \~german
 * @brief Programmstart
 *
 * \~english
 * @brief program start
 *
 * \~
 * @author Sebastian Pliet
 *
 */

#include "aetherise.h"
#include "cmd_line.h"
#include "utils.h"

#ifdef AETHER_WINDOWS
#include "windows.h"
#endif

#include <iostream>
#include <cstdlib>
#include <clocale>


using namespace aether;

const Version version {0,9,0};

const char* const ROOT_VERSION = "6.20/04";


void init_locale()
{
	char* locale = std::setlocale(LC_NUMERIC,""); // query user locale
	if (locale) {
		locale_german = starts_with(locale,"de") || starts_with(locale,"German");
	}
	std::setlocale(LC_NUMERIC,"C"); // back to normal
}



void show_german_help()
{
	std::cout << "aetherise [Schalter...] Dateien...\n\n";
	std::cout << "Ein Analysewerkzeug für die Datenblätter der Experimente von Dayton C. Miller\n";
	std::cout << "auf dem Mount Wilson in den Jahren 1925-1926.\n";
	std::cout << "Version " << version << "\n\n";
	std::cout << "Schalter:\n";
	std::cout << "  Filter werden mit Intervall [von,bis] angegeben.\n";
	std::cout << "  Ein fehlender Wert zeigt ein unbeschränktes Intervall an, z.B. [5.0,]\n";
	std::cout << "\n";
	std::cout << "  -no          <Intervall>  Nummer des Datenblattes\n";
	std::cout << "  -year        <Intervall>  Jahr des Datums\n";
	std::cout << "  -month       <Intervall>  Monat des Datums\n";
	std::cout << "  -day         <Intervall>  Tag des Datums\n";
	std::cout << "  -weight      <Intervall>  Gewicht\n";
	std::cout << "  -fringes     <Intervall>  Sichtbare Interferenzstreifen\n";
	std::cout << "  -time        <Intervall>  mittlere Beobachtungszeit\n";
	std::cout << "  -sidereal    <Intervall>  Sternzeit\n";
	std::cout << "  -T           <Intervall>  Temperatur\n";
	std::cout << "  -dT          <Intervall>  Temperaturänderung\n";
	std::cout << "  -mean_dT     <Intervall>  mittlere Temperaturänderung\n";
	std::cout << "  -TD          <Intervall>  maximale Temperaturdifferenz zu einem Zeitpunkt\n";
	std::cout << "  -adjust      <Intervall>  Anzahl der Neujustierungen\n";
    std::cout << "  -sign_correct             Vermerk 'sign correct' vorhanden\n";
    std::cout << "  -sign_correct_missing     Vermerk 'sign correct' fehlt\n";
	std::cout << "  -nw                       Schreibtisch in NW Ecke\n";
	std::cout << "  -sw                       Schreibtisch in SW Ecke\n";
	std::cout << "\n";
	std::cout << "  -amplitude   <Intervall>  Amplitude des Signals (Daten)\n";
	std::cout << "  -theory_amp  <Intervall>  Amplitude des Signals (Theorie)\n";
	std::cout << "  -drift       <Intervall>  Berichtigte mittlere Drift nach einer Umdrehung\n";
	std::cout << "  -uncertainty <Intervall>  Mittlere Standardunsicherheit\n";	
	std::cout << "\n";
	std::cout << "  Aktionen um bestimmte Daten zu erzeugen\n";
	std::cout << "\n";
	std::cout << "  -header                Metadaten ausgeben\n";
	std::cout << "  -raw                   Rohdaten\n";
	std::cout << "  -raw_reduced           Datenreduktion je Einzelmessung\n";
	std::cout << "  -reduce                Datenreduktion\n";	
	std::cout << "  -test                  Anderson-Darling-Test auf Normalverteilung\n";
	std::cout << "  -aggregate <Methode>   Datenblätter aggregieren\n";
	std::cout << "              list       Übersicht aller wichtigen Daten/Statistiken\n";
	std::cout << "              test       Anderson-Darling-Test auf Normalverteilung\n";		
	std::cout << "              mean       Mittelwerte gewählter Aktion\n";
	std::cout << "              sidereal   Amplituden je Sternzeit\n";
	std::cout << "              diff_chi   Ähnlichkeit aufeinander folgender Datenblätter,\n";
	std::cout << "                         ausgedrückt mittels Chi-Quadrat-Test.\n";
	//std::cout << "              model_chi  Güte der Übereinstimmung von Modell und Daten\n";
	//std::cout << "              params     Parameter für das Modell ermitteln\n";
	std::cout << "              signals    Datenblätter zur Signalextraktion auswählen\n";	
	std::cout << "              fit        Ausgleichsrechnung. Theorieparameter ermitteln\n";	
	std::cout << "\n";
	std::cout << "  Schalter von Aktionen\n";
	std::cout << "\n";	
	std::cout << "  -reduction <Methode>   Datenverarbeitungmethode\n";
	std::cout << "              Miller     Dayton Millers Methode (Voreinstellung)\n";	
	std::cout << "              separate   Jede Umdrehung einzeln reduzieren, dann mitteln\n";
	std::cout << "  -theory <Name>         Theorie\n";
	std::cout << "           classic       Klassische Äthertheorie\n";
	std::cout << "           ether         Äthertheorie von Lorentz (Voreinstellung)\n";
	std::cout << "           relativity    Spezielle Relativitätstheorie\n";
	std::cout << "  -single                Doppelperiode zu einer Periode mitteln\n";	
	std::cout << "  -subtract_theory       Theorie von Daten und Modell abziehen\n";
	std::cout << "  -add_theory            Theorie auf Daten addieren\n";
	std::cout << "  -invert_data           Vorzeichen der Daten ändern\n";
	std::cout << "  -invert_theory         Vorzeichen der Theorie ändern\n";
	std::cout << "  -invert_model          Vorzeichen des Fehlermodells ändern\n";
	std::cout << "  -invert_temp           Vorzeichen des Temperaturmodells ändern\n";
	std::cout << "  -data <Dateiname>      Daten im CSV-Format laden\n";
	std::cout << "  -subtract_data         Daten aus Datei von Messdaten abziehen.\n";
	std::cout << "  -subtract_model        Fehlermodell von Daten abziehen.\n";		
	std::cout << "  -disable_desk          Schreibtischsignal im Modell abschalten\n";
	std::cout << "  -disable_temp          Temperatursignal im Modell abschalten\n";
	std::cout << "  -disable_earth         Geschwindigkeitsvektor der Erdbahn nicht einrechnen\n";
	//std::cout << "  -epoch_params          Für jede Epoche und Schreibtischort einen\n";
	//std::cout << "                         gemeinsamen Satz an Modellparametern verwenden.\n";
	std::cout << "  -signals_dTD           Maximales ∆TD zwischen Datenblättern in °C\n";
	std::cout << "  -signals_ddT           Maximales ∆dT zwischen Datenblättern in °C\n";
	std::cout << "  -signals_dt            Minimales ∆t zwischen Datenblätteraggregaten in h\n";
	std::cout << "  -day_and_night         Signale auch aus Tag und Nacht Daten extrahieren\n";
	std::cout << "  -low_sun               Signale auch bei Sonnen-aufgang/untergang extrahieren\n";
	std::cout << "  -fit_amplitude         Bei der Ausgleichsrechnung (-aggregate fit)\n";
	std::cout << "                         nur die Amplitude nutzen, nicht das ganze Signal.\n";
	std::cout << "  -fit_sine              Anpassung an Phase/Amplitude von angepasstem Sinus\n";	
	std::cout << "  -fit_disable <Nummern> Signale abschalten. Nummern mit -stats sichtbar.\n";
	std::cout << "                         Erwartet Liste von Nummern durch Komma getrennt.\n";
	std::cout << "  -fix_ad                Parameter (α,δ) bei der Minimierung festhalten\n";
	std::cout << "  -minimizer <Name>      Minimierer für die Ausgleichsrechnung\n";
	std::cout << "              grad       Einfaches Gradientenverfahren (Voreinstellung)\n";
#ifdef AETHER_MINUIT
	std::cout << "              Minuit2    Minuit2 aus ROOT " << ROOT_VERSION << " vom CERN\n";
#endif
	std::cout << "  -delta_chi2 <Wert>     ∆χ² wird zur Berechnung der Parameterunsicherheiten\n";
	std::cout << "                         bei der Minimierung benutzt. Voreinstellung: 1.0\n";
	std::cout << "  -chi2_scale <Wert>     χ²-Skalierungsfaktor\n";
	std::cout << "  -theory_params <v,α,δ> Voreinstellung: (369000,11.2,-6.9) (m/s,h,°)\n";
	std::cout << "  -start_params <v,α,δ>  Startparameter für den Minimierer\n";
	std::cout << "  -n <Wert>              Brechungsindex auf festen Wert setzen\n";
	std::cout << "  -contour               Isolinien berechnen nach der Minimierung\n";
	std::cout << "  -ignore <Kürzel>       Angegebene Attribute ignorieren\n";
	std::cout << "   Eine Kombination von:\n";
	std::cout << "           -             Vorzeichen des Datenblattes umdrehen\n";
	std::cout << "           i             Vorzeichen einer Messung umdrehen\n";
	std::cout << "           b             Messung weglassen\n";	
	std::cout << "           r             Vorzeichen einer Messung umdrehen (Miller)\n";
	std::cout << "           c             Messung weglassen (Miller)\n";
	std::cout << "           R             Aufgehobenes r\n";
	std::cout << "           C             Aufgehobenes c\n";
	std::cout << "   Oder eine der Mengen:\n";
	std::cout << "           all           Alle Kürzel, außer r und c\n";
	std::cout << "           all!          Alle Kürzel\n";
	std::cout << "\n";
	std::cout << "  Andere Schalter\n";
	std::cout << "\n";
	std::cout << "  -validate              Widerspruchsfreiheit prüfen\n";
	std::cout << "\n";
	std::cout << "  Ausgabe\n";
	std::cout << "\n";
	std::cout << "  -stats                 Statistiken\n";
	std::cout << "  -model                 Modell anzeigen\n";
	std::cout << "  -no_data               Daten nicht ausgeben\n";
	std::cout << "  -no_theory             Theorie nicht ausgeben\n";
	std::cout << "  -csv                   Ausgabe im CSV-Format\n";	
	std::cout << "\n";
}



void show_english_help()
{
	std::cout << "aetherise [options...] files...\n\n";
	std::cout << "A tool to analyse the data sheets from the experiments of Dayton C. Miller\n";
	std::cout << "on top of Mount Wilson in the years 1925-1926.\n";
	std::cout << "Version " << version << "\n\n";
	std::cout << "Options:\n";
	std::cout << "  Filters have an interval [from,to].\n";
	std::cout << "  A missing value declares an unbounded interval, e.g. [5.0,]\n";
	std::cout << "\n";
	std::cout << "  -no          <interval>  Number of the data sheet\n";
	std::cout << "  -year        <interval>  Year of the date\n";
	std::cout << "  -month       <interval>  Month of the date\n";
	std::cout << "  -day         <interval>  Day of the date\n";
	std::cout << "  -weight      <interval>  Weight\n";
	std::cout << "  -fringes     <interval>  Visible fringes\n";
	std::cout << "  -time        <interval>  Time of observations (mean)\n";
	std::cout << "  -sidereal    <interval>  Sidereal time\n";
	std::cout << "  -T           <interval>  Temperature\n";
	std::cout << "  -dT          <interval>  Temperature changes\n";
	std::cout << "  -mean_dT     <interval>  Mean temperature change\n";
	std::cout << "  -TD          <interval>  Maximum temperature difference at a time\n";
	std::cout << "  -adjust      <interval>  Number of adjustments\n";
    std::cout << "  -sign_correct            Note 'sign correct' present\n";
    std::cout << "  -sign_correct_missing    Note 'sign correct' missing\n";
	std::cout << "  -nw                      Desk in NW corner\n";
	std::cout << "  -sw                      Desk in SW corner\n";
	std::cout << "\n";
	std::cout << "  -amplitude   <interval>  Amplitude of signal (reduced data)\n";
	std::cout << "  -theory_amp  <interval>  Amplitude of signal (theory)\n";
	std::cout << "  -drift       <interval>  Corrected mean drift after a turn\n";
	std::cout << "  -uncertainty <interval>  Mean standard uncertainty\n";	
	std::cout << "\n";
	std::cout << "  Actions to generate specific data\n";
	std::cout << "\n";		
	std::cout << "  -header                Show meta data\n";
	std::cout << "  -raw                   Raw data\n";
	std::cout << "  -raw_reduced           Data reduction for each turn\n";
	std::cout << "  -reduce                Data reduction\n";	
	std::cout << "  -test                  Anderson-Darling test for normality\n";
	std::cout << "  -aggregate <method>    Aggregate data sheets\n";
	std::cout << "              list       Overview of all relevant data\n";
	std::cout << "              test       Anderson-Darling test for normality\n";		
	std::cout << "              mean       Mean values of given action\n";
	std::cout << "              sidereal   Amplitudes at sidereal time\n";	
	std::cout << "              diff_chi   Similarity of sequenced data sheets,\n";
	std::cout << "                         using the chi squared test.\n";
	//std::cout << "              model_chi  Goodness of fit for the model\n";
	//std::cout << "              params     Find model parameters\n";
	std::cout << "              signals    Select data sheets for signal extraction\n";	
	std::cout << "              fit        Curve fitting. Find theory parameters\n";	
	std::cout << "\n";
	std::cout << "  Switches for actions\n";
	std::cout << "\n";	
	std::cout << "  -reduction <method>    Data reduction method\n";
	std::cout << "              Miller     Dayton Millers method (default)\n";	
	std::cout << "              separate   Reduce each turn separately, then average\n";
	std::cout << "  -theory <name>         Theory\n";
	std::cout << "           classic       Classic ether theory\n";
	std::cout << "           ether         Lorentz Ether Theory (default)\n";
	std::cout << "           relativity    Special Theory Of Relativity\n";
	std::cout << "  -single                Reduce double period to single period\n";	
	std::cout << "  -subtract_theory       Subtract theory from data and model\n";
	std::cout << "  -add_theory            Add theory to the data\n";
	std::cout << "  -invert_data           Changes the sign of the data\n";
	std::cout << "  -invert_theory         Changes the sign of the theory\n";
	std::cout << "  -invert_model          Changes the sign of the error model\n";
	std::cout << "  -invert_temp           Changes the sign of the temperature model\n";
	std::cout << "  -data <filename>       Load data in CSV format\n";
	std::cout << "  -subtract_data         Subtract data file from the data\n";
	std::cout << "  -subtract_model        Subtract error model from the data\n";		
	std::cout << "  -disable_desk          Disable desk signal in the model\n";
	std::cout << "  -disable_temp          Disable temperature signal in the model\n";	
	std::cout << "  -disable_earth         Do not add the earth orbital velocity vector\n";
	//std::cout << "  -epoch_params          Use one set of model parameters for each epoch\n";
	//std::cout << "                         and desk location\n";
	std::cout << "  -signals_dTD           Maximum ∆TD between data sheets in °C\n";
	std::cout << "  -signals_ddT           Maximum ∆dT between data sheets in °C\n";
	std::cout << "  -signals_dt            Minimum ∆t between data sheet aggregates in h\n";
	std::cout << "  -day_and_night         Extract signals from day and night data\n";
	std::cout << "  -low_sun               Extract signals at sunrise or sunset\n";
	std::cout << "  -fit_amplitude         Use only the amplitude, not the whole signal, \n";
	std::cout << "                         at minimization (-aggregate fit)\n";
	std::cout << "  -fit_sine              Fit against fitted sine curve\n";
	std::cout << "  -fit_disable <Numbers> Disable signals. Use -stats to find the numbers.\n";
	std::cout << "                         Expects a comma separated list of numbers.\n";
	std::cout << "  -fix_ad                Fix parameters (α,δ) at minimization\n";
	std::cout << "  -minimizer <name>      minimizer to use at minimization\n";
	std::cout << "              grad       Simple gradient descent (default)\n";
#ifdef AETHER_MINUIT
	std::cout << "              Minuit2    Minuit2 of ROOT " << ROOT_VERSION << " from CERN\n";
#endif
	std::cout << "  -delta_chi2 <value>    ∆χ² used for calculation of the uncertainties\n";
	std::cout << "                         of parameters at minimization. Default: 1.0\n";
	std::cout << "  -chi2_scale <value>    χ² scale factor\n";
	std::cout << "  -theory_params <v,α,δ> Default: (369000,11.2,-6.9) (m/s,h,°)\n";
	std::cout << "  -start_params <v,α,δ>  Start parameters for minimizing\n";
	std::cout << "  -n <value>             Set index of refraction to a fixed value\n";
	std::cout << "  -contour               Calculate contour after minimizing\n";
	std::cout << "  -ignore <code>         Ignore given attributes\n";
	std::cout << "   A combination of:\n";
	std::cout << "           -             Reverse sign of data sheet\n";
	std::cout << "           i             Change sign of a turn\n";
	std::cout << "           b             Cancel turn\n";	
	std::cout << "           r             Change sign of a turn (Miller)\n";
	std::cout << "           c             Cancel turn (Miller)\n";
	std::cout << "           R             Disabled r\n";
	std::cout << "           C             Disabled c\n";
	std::cout << "   Or one of the sets:\n";
	std::cout << "           all           All codes, without r and c\n";
	std::cout << "           all!          All codes\n";
	std::cout << "\n";
	std::cout << "  Other switches\n";
	std::cout << "\n";
	std::cout << "  -validate              Validate consistency\n";
	std::cout << "\n";
	std::cout << "  Output\n";
	std::cout << "\n";
	std::cout << "  -stats                 Stats\n";
	std::cout << "  -model                 Show model\n";
	std::cout << "  -no_data               Do not output data\n";
	std::cout << "  -no_theory             Do not output theory\n";
	std::cout << "  -csv                   Output in CSV format\n";	
	std::cout << "\n";
}


void show_help()
{
	if (locale_german)
		show_german_help();
	else
		show_english_help();
}



int main(int argc, char *argv[])
{
	init_locale();		

#ifdef AETHER_WINDOWS
	SetConsoleOutputCP(CP_UTF8);
	SetConsoleCP(CP_UTF8);
#endif

	if (argc <= 1) {
		show_help();
		return EXIT_FAILURE;
	}


	Filter filter;
	Action action = Action::Filename;
	bool aggregate = false;
	Options options;
	std::vector<std::string> filenames;

	parse_arguments(argc,argv,filter,action,aggregate,options,filenames);
	if (filenames.empty()) {
		std::cerr << "no input files\n";
		return EXIT_FAILURE;
	}

	std::vector<DataSheet> data_sheets;
	int n = 0;

	try {
		for (const std::string& filename : filenames) {		
			DataSheet data_sheet = load_data_sheet_csv(filename,options);
			if (selected(data_sheet,options,filter)) {
				n++;
				if (options.validate)
					validate(data_sheet,filename,std::cerr);
	
				if (aggregate)
					data_sheets.push_back(std::move(data_sheet));
				else
					execute(action,options,data_sheet,filename);
			}
		}
	
		if (n == 0) {
			std::cerr << "no data sheets selected\n";
			return EXIT_FAILURE;
		}
	
		if (aggregate)
			execute(action,options,data_sheets);			
	}
	catch(ExitException& e) {
		if (!e.message.empty())
			std::cerr << e.message << "\n";
		std::exit(EXIT_FAILURE);
	}

	return EXIT_SUCCESS;
}

