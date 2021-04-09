#include "models.h"

#include "mathematics.h"
#include "data_reduction.h"
#include "physics.h"

#include <algorithm>
#include <iostream>
#include <sstream>

namespace aether {





std::array<double,16> T_distribution(const DataSheet::Thermometers& thermometers)
{
	std::array<double,16> Ts;

	for (int i=1; i<17; i++) { // azimuth
		if (i>=1 && i<5) {
			double d = thermometers.E - thermometers.N;
			Ts.at(i-1) = thermometers.N + d/4*(i-1);
		}
		else if (i>=5 && i<9) {
			double d = thermometers.S - thermometers.E;
			Ts.at(i-1) = thermometers.E + d/4*(i-5);
		}
		else if (i>=9 && i<13) {
			double d = thermometers.W - thermometers.S;
			Ts.at(i-1) = thermometers.S + d/4*(i-9);
		}
		else {
			double d = thermometers.N - thermometers.W;
			Ts.at(i-1) = thermometers.W + d/4*(i-13);
		}
	}

	// some non linearity
	Ts.at(1) = (Ts.at(0)+Ts.at(1))/2;
	Ts.at(15) = (Ts.at(0)+Ts.at(15))/2;
	Ts.at(3) = (Ts.at(4)+Ts.at(3))/2;
	Ts.at(5) = (Ts.at(4)+Ts.at(5))/2;
	Ts.at(7) = (Ts.at(8)+Ts.at(7))/2;
	Ts.at(9) = (Ts.at(8)+Ts.at(9))/2;
	Ts.at(11) = (Ts.at(12)+Ts.at(11))/2;
	Ts.at(13) = (Ts.at(12)+Ts.at(13))/2;

	return Ts;
}





optional<DataSheet::Thermometers> in_Kelvin(const optional<DataSheet::Thermometers>& th)
{
	if (!th.has_value())
		return th;

	DataSheet::Thermometers k;

	k.N = in_Kelvin(th->N);
	k.E = in_Kelvin(th->E);
	k.S = in_Kelvin(th->S);
	k.W = in_Kelvin(th->W);

	return {k};
}


void warm_up(double& Tx, double T,double )
{	
	Tx += (T-Tx)/2;	
}



std::array<double,17>
systematic_error_displacements(const optional<DataSheet::Thermometers>& start,
							   const optional<DataSheet::Thermometers>& end,
							   const ErrorModelParameters& params,
							   const Options& options,bool desk_disabled)
{
	const double pi = AETHER_PI;

	std::array<double,16> Ds {};
	const double desk = params.desk;
	const int desk_azimuth = std::floor(params.desk_azimuth);	

	double az = (params.desk_azimuth-desk_azimuth) * pi/8;
	// create non linear desk effect distribution
	Ds.at((desk_azimuth-1+15)%16) += desk * std::sin(pi/2-pi/4-az); // +15=-1
	Ds.at((desk_azimuth-1   )%16) += desk * std::sin(pi/2     -az);
	Ds.at((desk_azimuth-1+1 )%16) += desk * std::sin(pi/2+pi/4-az);
	Ds.at((desk_azimuth-1+2 )%16) += desk * std::sin(pi/2+pi/2-az);


	std::array<double,17> dts {};
	// desk effect
	if (options.enable_desk && !desk_disabled) {
		for (int i=0;i<17;i++) {
			double D0 = Ds.at((i+8 )%16); // S (optic)
			double D1 = Ds.at((i   )%16); // N (main mirror)
			double D2 = Ds.at((i+4 )%16); // E
			double D3 = Ds.at((i+12)%16); // W

			double t1 = D0 + D1;
			double t2 = D2 + D3;
			double dt = t2-t1;

			dts.at(i) += -0.5*dt*(options.invert_model ? -1 : 1); // 0.5 because of old params for fringes, not for wave length
		}
	}

	if (options.enable_temp) {
		if (start.has_value() || end.has_value()) {

			const auto meanth = *mean_thermometers(in_Kelvin(start),in_Kelvin(end));  // value always exists
			const auto meanT = mean_T(meanth);
			std::array<double,16> Ts = T_distribution(meanth);
			std::rotate(Ts.begin(),Ts.begin()+1,Ts.end()); // house is rotated 20° from N

			double Tt0,Tt1,Tt2,Tt3;
			Tt0=Tt1=Tt2=Tt3=sqr(meanT);

			// warm up
			int n=1;
			for (int k=0;k<n;k++) {
				for (int i=0;i<16;i++) {
					double T0 = Ts.at((i+8 )%16); // S (optic)
					double T1 = Ts.at((i   )%16); // N (main mirror)
					double T2 = Ts.at((i+4 )%16); // E
					double T3 = Ts.at((i+12)%16); // W

					// increase differences non linear
					T0=sqr(T0);
					T1=sqr(T1);
					T2=sqr(T2);
					T3=sqr(T3);

					warm_up(Tt0,T0,0);
					warm_up(Tt1,T1,0);
					warm_up(Tt2,T2,0);
					warm_up(Tt3,T3,0);
				}
			}


			// temperature effect
			for (int i=0;i<17;i++) {
				double T0 = Ts.at((i+8 )%16); // S (optic)
				double T1 = Ts.at((i   )%16); // N (main mirror)
				double T2 = Ts.at((i+4 )%16); // E
				double T3 = Ts.at((i+12)%16); // W

				// increase differences non linear
				T0=sqr(T0);
				T1=sqr(T1);
				T2=sqr(T2);
				T3=sqr(T3);

				warm_up(Tt0,T0,0);
				warm_up(Tt1,T1,0);
				warm_up(Tt2,T2,0);
				warm_up(Tt3,T3,0);

				double t1 = (Tt0-meanT) + (Tt1-meanT);
				double t2 = (Tt2-meanT) + (Tt3-meanT);
				double dt = (t2-t1)*0.0025;

				dts.at(i) += -0.5*dt*params.temperature*(options.invert_model!=options.invert_temp ? -1 : 1);
			}
		}		
	}

	if (options.single)
		reduce_to_single_period(dts);

	// needed if temperature model is not symmetric
	double mean_ordinate = std::accumulate(dts.begin(),dts.end()-1,0.0) / 16.0;
	for (double& d : dts) {
		d -= mean_ordinate;
	}

	return dts;
}


bool in_group(const DataSheet& data_sheet,int from, int to)
{
	return data_sheet.no>=from && data_sheet.no<=to;
}


// Used options: -single

// Abbreviations
//-------------------------------------------------------------
// Weather (W): c=cloud, f=fog, w=wind, r=rain, (empty)=clear
// Temperature Difference (TD): s=significant signal


// Group		W	t		T		TD		dT			check
//------------------------------------------------
// Apr 103-106		16:30	13		1.7	s	-			amplitude-
// Apr 107-108	f	15:30	4		0.4		-0.5		!
// Apr 110-111	f	 0:40	0.6		0.4		 0.5		! (desk effect gone?)
// Apr 113-116		 3:00	1.3		0.5	s				!
// Apr 122-125	f	16:10	10		0.5	s	-0.4		amplitude-
// Apr 127-128		16:00	16		1.6	s	-0.2		!
// Apr 130-133		 1:30	9		0.8	s	 0.2		! (desk effect gone?)
//------------------------------------------------
//const ErrorModelParameters params_apr_nw_103 {14.75, 0.065, 0.045};
const ErrorModelParameters params_apr_nw_103 {12.75, 0.195, 0.08};
//const ErrorModelParameters params_apr_nw_107 {14.25, 0.025, 0.02}; // (r)
const ErrorModelParameters params_apr_nw_107 {12.75, 0.06, 0.0};
const ErrorModelParameters params_apr_nw_110 {14.5, 0.01, 0.02};// anomaly
//const ErrorModelParameters params_apr_nw_113 {14.5, 0.045, 0.025};// (r)
const ErrorModelParameters params_apr_nw_113 {12, 0.065, 0.0};
//const ErrorModelParameters params_apr_nw_122 {14.5, 0.04, 0.03};// (r)
//const ErrorModelParameters params_apr_nw_122 {12.75, 0.045, -0.035};
const ErrorModelParameters params_apr_nw_122 {12.75, 0.1, 0.04};
//const ErrorModelParameters params_apr_nw_127 {14.5, 0.03, 0.03}; // (r)
//const ErrorModelParameters params_apr_nw_127 {13.25, 0.11, 0.0};
const ErrorModelParameters params_apr_nw_127 {13.5, 0.045, -0.06};
const ErrorModelParameters params_apr_nw_130 {14.5, 0.01, 0.02};// anomaly
const ErrorModelParameters params_apr_nw {14.5, 0.027, 0.02};




// Group		W	t		T		TD		dT			check
//------------------------------------------------
// Apr 117-120		5:00	1.3		0.5		0			!
//------------------------------------------------
// 119-120 sunrise, temperature invalid?
//const ErrorModelParameters params_apr_sw_117 {11.25, 0.135, 0.01};
const ErrorModelParameters params_apr_sw_117 {11.75, 0.045, -0.02};
const ErrorModelParameters params_apr_sw {11.25, 0.135, 0.01};



// Group		W	t		T		TD		dT			check
//------------------------------------------------
// Aug  4-10		 3:30	20		0.6		-0.2		!
// Aug 12-17		 3:30	19		0.4		-0.2		!
// Aug 21-26		 3:30	19		0.5		-0.2		?
// Aug 29-30		 8:00	27		0.7		 0.6		phase-1
// Aug 31-36		16:00	28		0.4		-0.4		!
// Aug 37-40		21:00	23		0.3		-0.5		!
// Aug 45-48		17:00	28		1.2	s	-0.4		?
// Aug 49-54		16:00	29		0.5		-0.2		!
// Aug 60-63		14:00	26		0.5		 0.1		phase+1
//------------------------------------------------
const ErrorModelParameters params_aug_nw_4 {14, 0.03, 0.025};
const ErrorModelParameters params_aug_nw_12 {13.75, 0.02, 0.025};
const ErrorModelParameters params_aug_nw_21 {12, 0.145, 0.005};
//const ErrorModelParameters params_aug_nw_21 {16, 0.035, 0.03};
const ErrorModelParameters params_aug_nw_29 {13.25, 0.045, 0.035};
const ErrorModelParameters params_aug_nw_31 {12.75, 0.045, 0.035};// anomaly (without 32-34, distorted?)
const ErrorModelParameters params_aug_nw_37 {12.75, 0.04, 0.045};
const ErrorModelParameters params_aug_nw_45 {12.75, 0.075, 0.01};
const ErrorModelParameters params_aug_nw_49 {13, 0.04, 0.045}; // anomaly (without 50-52, distorted?)
const ErrorModelParameters params_aug_nw_60 {12, 0.03, 0.03};
const ErrorModelParameters params_aug_nw {12.75, 0.075, 0.01};



// Group		W	t		T		TD		dT			check
//------------------------------------------------
// Aug 69-73		11:21	23		0.7 s	 0.4		!
// Aug 74-77		14:30	26		0.4		 0.1		phase-1 (76-77)
// Aug 80-85		 7:00	18		0.7	s	 0.4		!
// Aug 86-87		10:00	25		0.5		 0.2		!
// Aug 88-90		 0:00	17		0.4	s	-0.1		!
// Aug 91-96		 7:00	17		0.8	s	 0.5		!
//------------------------------------------------
const ErrorModelParameters params_aug_sw_69 {9, 0.055, 0.045};
const ErrorModelParameters params_aug_sw_74 {10, 0.06, 0.05};
const ErrorModelParameters params_aug_sw_78 {11, 0.04, 0.03};
const ErrorModelParameters params_aug_sw_80 {9.75, 0.075, 0.035};
const ErrorModelParameters params_aug_sw_86 {9.25, 0.05, 0.055};
const ErrorModelParameters params_aug_sw_88 {9.25, 0.045, 0.075};
const ErrorModelParameters params_aug_sw_91 {9.5, 0.06, 0.05};
const ErrorModelParameters params_aug_sw {9.75, 0.075, 0.035};



// Group		W	t		T		TD		dT			check
//------------------------------------------------
// Sep 1-3			17:00	23		0.2		-0.2		!
// Sep 14-19	f	15:30	17		0.4	s	±0.1		anomaly
// Sep 23-25		17:00	22		0.2		-0.4		!
// Sep 26-28		14:30	26		0.1		 0.2		!
// Sep 39-41		 4:30	13		0.5		-0.2		!
// Sep 49-56	c	22:30	15		0.4		-0			!
// Sep 57-62	c	10:30	14		0.5	s	 0.6		!
// Sep 63-64		17:00	16		0.5		-0.6		!
// Sep 75-83		 4:00	14		0.3	s	-0			! (sw?)
//------------------------------------------------
const ErrorModelParameters params_sep_nw_1 {12.5, 0.035, 0.04};
const ErrorModelParameters params_sep_nw_8 {13, 0.035, 0.045};
const ErrorModelParameters params_sep_nw_14 {12.75, 0.045, 0.055};// anomaly
const ErrorModelParameters params_sep_nw_20_21 {12.75, 0.045, 0.035};
const ErrorModelParameters params_sep_nw_22 {12.25, 0.03, 0.035};
const ErrorModelParameters params_sep_nw_26 {12.25, 0.045, 0.035};
const ErrorModelParameters params_sep_nw_29 {13, 0.055, 0.045};
const ErrorModelParameters params_sep_nw_35 {12.75, 0.035, 0.03}; // 39-44
const ErrorModelParameters params_sep_nw_49 {13, 0.045, 0.02};
const ErrorModelParameters params_sep_nw_57 {12.75, 0.03, 0.03};
const ErrorModelParameters params_sep_nw_63 {13, 0.04, 0.03};
//const ErrorModelParameters params_sep_nw_75 {14.25, 0.05, 0.055}; // anomaly (r)
const ErrorModelParameters params_sep_nw {12.75, 0.045, 0.025};



// Group		W	t		T		TD		dT			check
//------------------------------------------------
// Sep 75-83		4:00	14		0.3	s	-0			!
//------------------------------------------------
// manually marked as SW, not by Miller
const ErrorModelParameters params_sep_sw_75 {10.75, 0.04, 0.05}; // 78-81
const ErrorModelParameters params_sep_sw = {11, 0.045, 0.035};



// Group		W	t		T		TD		dT			check
//------------------------------------------------
// Feb 18-20	c	17:30	13		0.4		-0.2		desk effect gone?
// Feb 43-46	c	18:00	13		0.3		-0.3		desk effect gone?
// Feb 48-52		23:00	11		0.1		 0.1		!
// Feb 53-58		 6:00	 9		0.4	s	 0.1		!
// Feb 69-73		 8:00	 9		0.5	s	 0.2		!
//------------------------------------------------
const ErrorModelParameters params_feb_nw_18 {14, 0, 0};// 18-20 (D)
//const ErrorModelParameters params_feb_nw_21 {13, 0.035, 0.045};//(r)
const ErrorModelParameters params_feb_nw_21 {15, 0.04, 0.07};
const ErrorModelParameters params_feb_nw_32 {13.25, 0.04, 0.05};//(r)
const ErrorModelParameters params_feb_nw_43 {14, 0, 0};// 43-46 (D)
const ErrorModelParameters params_feb_nw_47 {14.75, 0.045, 0};// 47-52 (r)
//const ErrorModelParameters params_feb_nw_53 {14.25, 0.04, 0.03};// 53-58 (r)
const ErrorModelParameters params_feb_nw_53 {14, 0.04, 0.02};// 53-58 (r)
const ErrorModelParameters params_feb_nw_69 {14, 0.04, 0.03};// 69-73 (r)
const ErrorModelParameters params_feb_nw {14.25, 0.04, 0.03};


// Group		W	t		T		TD		dT			check
//------------------------------------------------
// Feb 74-79	c	11:00	12		0.9		 0.2		!
// Feb 80-91		 3:40	6		0.6	s	-0.1		!
// Feb 92-94	f	16:00	7		0.3		-0			!
// Feb 95-101	wr	 2:00	4		0.2		-0			!
//------------------------------------------------
const ErrorModelParameters params_feb_sw_74 {8.5, 0.03, 0.03}; // 74-79
const ErrorModelParameters params_feb_sw_80 {10, 0.035, 0.04}; // 80-91
const ErrorModelParameters params_feb_sw_92 {10.25, 0.045, 0.035}; // 92-94
const ErrorModelParameters params_feb_sw_95 {9.75, 0.04, 0.045}; // 95-101 (without 97, 100)
const ErrorModelParameters params_feb_sw {10, 0.035, 0.04};





std::array<double,17>
systematic_error_displacements(const DataSheet& data_sheet,const Options& options)
{	
	ErrorModelParameters params;

	switch(data_sheet.date.month) { // epochs
	case 2:
		params = data_sheet.desk_in_sw ? params_feb_sw : params_feb_nw;

		if (!options.epoch_params) {
			if (in_group(data_sheet,18,20))
				params = params_feb_nw_18;
			else if (in_group(data_sheet,21,26))
				params = params_feb_nw_21;
			else if (in_group(data_sheet,32,36))
				params = params_feb_nw_32;
			else if (in_group(data_sheet,43,46))
				params = params_feb_nw_43;
			else if (in_group(data_sheet,47,52))
				params = params_feb_nw_47;
			else if (in_group(data_sheet,53,58))
				params = params_feb_nw_53;
			else if (in_group(data_sheet,69,73))
				params = params_feb_nw_69;
			else if (in_group(data_sheet,74,79))
				params = params_feb_sw_74;
			else if (in_group(data_sheet,80,91))
				params = params_feb_sw_80;
			else if (in_group(data_sheet,92,94))
				params = params_feb_sw_92;
			else if (in_group(data_sheet,95,101))
				params = params_feb_sw_95;
			else {
				std::cerr << "WARNING: no model parameters for sheet " << data_sheet.no
						  << " found, using default\n";
			}
		}
		break;
	case 3:
	case 4:
		params = data_sheet.desk_in_sw ? params_apr_sw : params_apr_nw;

		if (!options.epoch_params) {
			if (in_group(data_sheet,103,106))
				params = params_apr_nw_103;
			else if (in_group(data_sheet,107,108))
				params = params_apr_nw_107;
			else if (in_group(data_sheet,110,111))
				params = params_apr_nw_110;
			else if (in_group(data_sheet,113,116))
				params = params_apr_nw_113;
			else if (in_group(data_sheet,117,120))
				params = params_apr_sw_117;
			else if (in_group(data_sheet,122,126))
				params = params_apr_nw_122;
			else if (in_group(data_sheet,127,129))
				params = params_apr_nw_127;
			else if (in_group(data_sheet,130,133))
				params = params_apr_nw_130;
			else {
				std::cerr << "WARNING: no model parameters for sheet " << data_sheet.no
						  << " found, using default\n";
			}
		}
		break;
	case 7:
	case 8:
		params = data_sheet.desk_in_sw ? params_aug_sw : params_aug_nw;

		if (!options.epoch_params) {
			if (in_group(data_sheet,4,10))
				params = params_aug_nw_4;
			else if (in_group(data_sheet,12,19))
				params = params_aug_nw_12;
			else if (in_group(data_sheet,21,28))
				params = params_aug_nw_21;
			else if (in_group(data_sheet,29,30))
				params = params_aug_nw_29;
			else if (in_group(data_sheet,31,36))
				params = params_aug_nw_31;
			else if (in_group(data_sheet,37,40))
				params = params_aug_nw_37;
			else if (in_group(data_sheet,45,48))
				params = params_aug_nw_45;
			else if (in_group(data_sheet,49,54))
				params = params_aug_nw_49;
			else if (in_group(data_sheet,60,63))
				params = params_aug_nw_60;
			else if (in_group(data_sheet,69,73))
				params = params_aug_sw_69;
			else if (in_group(data_sheet,74,77))
				params = params_aug_sw_74;
			else if (in_group(data_sheet,78,79))
				params = params_aug_sw_78;
			else if (in_group(data_sheet,80,85))
				params = params_aug_sw_80;
			else if (in_group(data_sheet,86,87))
				params = params_aug_sw_86;
			else if (in_group(data_sheet,88,90))
				params = params_aug_sw_88;
			else if (in_group(data_sheet,91,96))
				params = params_aug_sw_91;
			else {
				std::cerr << "WARNING: no model parameters for sheet " << data_sheet.no
						  << " found, using default\n";
			}
		}
		break;
	case 9:
		params = data_sheet.desk_in_sw ? params_sep_sw : params_sep_nw;

		if (!options.epoch_params) {
			if (in_group(data_sheet,1,3))
				params = params_sep_nw_1;
			else if (in_group(data_sheet,8,13))
				params = params_sep_nw_8;
			else if (in_group(data_sheet,14,19))
				params = params_sep_nw_14;
			else if (in_group(data_sheet,20,21))
				params = params_sep_nw_20_21;
			else if (in_group(data_sheet,22,25))
				params = params_sep_nw_22;
			else if (in_group(data_sheet,26,28))
				params = params_sep_nw_26;
			else if (in_group(data_sheet,29,34))
				params = params_sep_nw_29;
			else if (in_group(data_sheet,35,44))
				params = params_sep_nw_35;
			else if (in_group(data_sheet,49,56))
				params = params_sep_nw_49;
			else if (in_group(data_sheet,57,62))
				params = params_sep_nw_57;
			else if (in_group(data_sheet,63,65))
				params = params_sep_nw_63;
			else if (in_group(data_sheet,75,83))
				params = params_sep_sw_75;
			else {
				std::cerr << "WARNING: no model parameters for sheet " << data_sheet.no
						  << " found, using default\n";
			}
		}
		break;
	default:		
		throw std::runtime_error(SSSTR("unknown epoch (month) " << data_sheet.date.month));
	}

	if (options.enable_temp) {
		if (!data_sheet.thermometers_start.has_value() && !data_sheet.thermometers_end.has_value()) {
			std::cerr << "WARNING: no temperature model for sheet " << data_sheet.no << "\n";
		}
	}

	return systematic_error_displacements(data_sheet.thermometers_start,data_sheet.thermometers_end,
										  params,options,data_sheet.desk_disabled);
}







}//
