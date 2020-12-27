
#include "test.h"

#include "../DataSheet.h"

#include <algorithm>

using namespace aether;



TEST(DataSheet,sort)
{

	{ // already sorted
		DataSheet ds1; ds1.date = {2000,1,1}; ds1.no = 10;
		DataSheet ds2; ds2.date = {2000,1,1}; ds2.no = 20;
		DataSheet ds3; ds3.date = {2000,1,2}; ds3.no = 0;
		DataSheet ds4; ds4.date = {2000,1,2}; ds4.no = 13;
		DataSheet ds5; ds5.date = {2000,1,2}; ds5.no = 20;

		std::vector<DataSheet> dss;
		dss.push_back(ds1);
		dss.push_back(ds2);
		dss.push_back(ds3);
		dss.push_back(ds4);
		dss.push_back(ds5);

		std::sort(dss.begin(),dss.end());

		ASSERT_TRUE(dss.at(0).date==ds1.date && dss.at(0).no==ds1.no);
		ASSERT_TRUE(dss.at(1).date==ds2.date && dss.at(1).no==ds2.no);
		ASSERT_TRUE(dss.at(2).date==ds3.date && dss.at(2).no==ds3.no);
		ASSERT_TRUE(dss.at(3).date==ds4.date && dss.at(3).no==ds4.no);
		ASSERT_TRUE(dss.at(4).date==ds5.date && dss.at(4).no==ds5.no);
	}


	{ // already sorted backwards
		DataSheet ds1; ds1.date = {2000,1,1}; ds1.no = 10;
		DataSheet ds2; ds2.date = {2000,1,1}; ds2.no = 20;
		DataSheet ds3; ds3.date = {2000,1,2}; ds3.no = 0;
		DataSheet ds4; ds4.date = {2000,1,2}; ds4.no = 13;
		DataSheet ds5; ds5.date = {2000,1,2}; ds5.no = 20;

		std::vector<DataSheet> dss;
		dss.push_back(ds5);
		dss.push_back(ds4);
		dss.push_back(ds3);
		dss.push_back(ds2);
		dss.push_back(ds1);

		std::sort(dss.begin(),dss.end());

		ASSERT_TRUE(dss.at(0).date==ds1.date && dss.at(0).no==ds1.no);
		ASSERT_TRUE(dss.at(1).date==ds2.date && dss.at(1).no==ds2.no);
		ASSERT_TRUE(dss.at(2).date==ds3.date && dss.at(2).no==ds3.no);
		ASSERT_TRUE(dss.at(3).date==ds4.date && dss.at(3).no==ds4.no);
		ASSERT_TRUE(dss.at(4).date==ds5.date && dss.at(4).no==ds5.no);
	}


	{ // unsorted
		DataSheet ds3; ds3.date = {2000,1,2}; ds3.no = 0;
		DataSheet ds2; ds2.date = {2000,1,1}; ds2.no = 20;
		DataSheet ds5; ds5.date = {2000,1,2}; ds5.no = 20;
		DataSheet ds4; ds4.date = {2000,1,2}; ds4.no = 13;
		DataSheet ds1; ds1.date = {2000,1,1}; ds1.no = 10;


		std::vector<DataSheet> dss;
		dss.push_back(ds3);
		dss.push_back(ds2);
		dss.push_back(ds5);
		dss.push_back(ds4);
		dss.push_back(ds1);

		std::sort(dss.begin(),dss.end());

		ASSERT_TRUE(dss.at(0).date==ds1.date && dss.at(0).no==ds1.no);
		ASSERT_TRUE(dss.at(1).date==ds2.date && dss.at(1).no==ds2.no);
		ASSERT_TRUE(dss.at(2).date==ds3.date && dss.at(2).no==ds3.no);
		ASSERT_TRUE(dss.at(3).date==ds4.date && dss.at(3).no==ds4.no);
		ASSERT_TRUE(dss.at(4).date==ds5.date && dss.at(4).no==ds5.no);
	}
}



TEST(DataSheet_Date,operator_equal)
{
	{
		DataSheet::Date a {1925,02,01};
		DataSheet::Date b {1925,02,01};

		ASSERT_TRUE(a==b);
	}

	{
		DataSheet::Date a {1925,02,01};
		DataSheet::Date b {1925,02,02};

		ASSERT_FALSE(a==b);
	}

	{
		DataSheet::Date a {1925,02,01};
		DataSheet::Date b {1925,01,01};

		ASSERT_FALSE(a==b);
	}

	{
		DataSheet::Date a {1925,02,01};
		DataSheet::Date b {1926,02,01};

		ASSERT_FALSE(a==b);
	}
}


TEST(DataSheet_Date,operator_less)
{
	{
		DataSheet::Date a {1925,02,01};
		DataSheet::Date b {1925,02,01};

		ASSERT_FALSE(a<b);
	}

	{
		DataSheet::Date a {1925,02,01};
		DataSheet::Date b {1925,02,02};

		ASSERT_TRUE(a<b);
	}

	{
		DataSheet::Date a {1925,02,02};
		DataSheet::Date b {1925,02,01};

		ASSERT_FALSE(a<b);
	}

	{
		DataSheet::Date a {1925,02,01};
		DataSheet::Date b {1925,03,01};

		ASSERT_TRUE(a<b);
	}

	{
		DataSheet::Date a {1925,02,02};
		DataSheet::Date b {1925,01,01};

		ASSERT_FALSE(a<b);
	}

	{
		DataSheet::Date a {1925,02,01};
		DataSheet::Date b {1926,02,01};

		ASSERT_TRUE(a<b);
	}

	{
		DataSheet::Date a {1926,02,01};
		DataSheet::Date b {1925,02,01};

		ASSERT_FALSE(a<b);
	}

	{
		DataSheet::Date a {1926,01,02};
		DataSheet::Date b {1925,02,01};

		ASSERT_FALSE(a<b);
	}
}




TEST(DataSheet_Time,operator_equal)
{
	{
		DataSheet::Time a {02,01};
		DataSheet::Time b {02,01};

		ASSERT_TRUE(a==b);
	}

	{
		DataSheet::Time a {02,01};
		DataSheet::Time b {02,02};

		ASSERT_FALSE(a==b);
	}

	{
		DataSheet::Time a {02,01};
		DataSheet::Time b {01,01};

		ASSERT_FALSE(a==b);
	}

	{
		DataSheet::Time a {02,01};
		DataSheet::Time b {01,02};

		ASSERT_FALSE(a==b);
	}
}





TEST(DataSheet_Time,operator_less)
{
	{
		DataSheet::Time a {02,01};
		DataSheet::Time b {02,01};

		ASSERT_FALSE(a<b);
	}


	{
		DataSheet::Time a {02,01};
		DataSheet::Time b {02,02};

		ASSERT_TRUE(a<b);
	}

	{
		DataSheet::Time a {02,02};
		DataSheet::Time b {02,01};

		ASSERT_FALSE(a<b);
	}


	{
		DataSheet::Time a {01,02};
		DataSheet::Time b {02,01};

		ASSERT_TRUE(a<b);
	}


	{
		DataSheet::Time a {02,01};
		DataSheet::Time b {01,02};

		ASSERT_FALSE(a<b);
	}
}



TEST(time_to_h,test)
{
	{
		ASSERT_TRUE(time_to_h({0,00}) == 0);
		ASSERT_TRUE(time_to_h({24,00}) == 24);
		ASSERT_TRUE(time_to_h({-24,00}) == -24);
		ASSERT_TRUE(time_to_h({12,00}) == 12);
		ASSERT_TRUE(time_to_h({12,30}) == 12.5);
		ASSERT_TRUE(time_to_h({12,59}) == 12+59.0/60.0);
	}
}

TEST(h_to_time,test)
{
	{
		DataSheet::Time t = h_to_time(0);
		DataSheet::Time x {0,00};
		ASSERT_TRUE(t==x);
	}

	{
		DataSheet::Time t = h_to_time(24);
		DataSheet::Time x {0,00};
		ASSERT_TRUE(t==x);
	}

	{
		DataSheet::Time t = h_to_time(-24);
		DataSheet::Time x {0,00};
		ASSERT_TRUE(t==x);
	}

	{
		DataSheet::Time t = h_to_time(12);
		DataSheet::Time x {12,00};
		ASSERT_TRUE(t==x);
	}

	{
		DataSheet::Time t = h_to_time(-12);
		DataSheet::Time x {-12,00};
		ASSERT_TRUE(t==x);
	}

	{
		DataSheet::Time t = h_to_time(12.5);
		DataSheet::Time x {12,30};
		ASSERT_TRUE(t==x);
	}

	{
		DataSheet::Time t = h_to_time(12.999);
		DataSheet::Time x {12,59};
		ASSERT_TRUE(t==x);
	}
}


TEST(max_T,test)
{
	{
		DataSheet::Thermometers th {0,0,0,0};
		auto m = max_T(th);
		ASSERT_TRUE(m==0);
	}

	{
		DataSheet::Thermometers th {1,3,4,2};
		auto m = max_T(th);
		ASSERT_TRUE(m==4);
	}

	{
		DataSheet::Thermometers th {1,3,-4,2};
		auto m = max_T(th);
		ASSERT_TRUE(m==3);
	}

	{
		DataSheet::Thermometers th {4,3,-4,2};
		auto m = max_T(th);
		ASSERT_TRUE(m==4);
	}

	{
		DataSheet::Thermometers th {-1,-3,-4,2};
		auto m = max_T(th);
		ASSERT_TRUE(m==2);
	}
}



TEST(min_T,test)
{
	{
		DataSheet::Thermometers th {0,0,0,0};
		auto m = min_T(th);
		ASSERT_TRUE(m==0);
	}

	{
		DataSheet::Thermometers th {1,2,3,4};
		auto m = min_T(th);
		ASSERT_TRUE(m==1);
	}

	{
		DataSheet::Thermometers th {2,1,4,3};
		auto m = min_T(th);
		ASSERT_TRUE(m==1);
	}

	{
		DataSheet::Thermometers th {2,4,1,3};
		auto m = min_T(th);
		ASSERT_TRUE(m==1);
	}

	{
		DataSheet::Thermometers th {4,3,2,1};
		auto m = min_T(th);
		ASSERT_TRUE(m==1);
	}

	{
		DataSheet::Thermometers th {-4,3,2,1};
		auto m = min_T(th);
		ASSERT_TRUE(m==-4);
	}
}



TEST(mean_T,test)
{
	{
		ASSERT_TRUE(mean_T({0,0,0,0}) == 0);
		ASSERT_TRUE(mean_T({-1,2,-2,1}) == 0);
		ASSERT_TRUE(mean_T({0.5,1.5,4,-2}) == 1);
		ASSERT_TRUE(mean_T({2,2,2,2}) == 2);
	}
}


TEST(mean_dT,test)
{
	{
		ASSERT_TRUE(mean_dT({1,1,1,1},{2,2,0,0}) == 0);
		ASSERT_TRUE(mean_dT({1,1,1,1},{3,3,1,1}) == 1);
		ASSERT_TRUE(mean_dT({3,3,1,1},{1,1,1,1}) == -1);
	}
}



TEST(time_of_day,datasheet)
{
	{
		DataSheet data_sheet;
		data_sheet.mean_observation_time = {7,00};
		data_sheet.date.month = 2;
		ASSERT_TRUE(time_of_day(data_sheet) == TimeOfDay::Sunrise);
	}

	{
		DataSheet data_sheet;
		data_sheet.mean_observation_time = {6,00};
		data_sheet.date.month = 3;
		ASSERT_TRUE(time_of_day(data_sheet) == TimeOfDay::Sunrise);
	}

	{
		DataSheet data_sheet;
		data_sheet.mean_observation_time = {6,00};
		data_sheet.date.month = 4;
		ASSERT_TRUE(time_of_day(data_sheet) == TimeOfDay::Sunrise);
	}

	{
		DataSheet data_sheet;
		data_sheet.mean_observation_time = {5,30};
		data_sheet.date.month = 7;
		ASSERT_TRUE(time_of_day(data_sheet) == TimeOfDay::Sunrise);
	}

	{
		DataSheet data_sheet;
		data_sheet.mean_observation_time = {5,30};
		data_sheet.date.month = 8;
		ASSERT_TRUE(time_of_day(data_sheet) == TimeOfDay::Sunrise);
	}

	{
		DataSheet data_sheet;
		data_sheet.mean_observation_time = {6,00};
		data_sheet.date.month = 9;
		ASSERT_TRUE(time_of_day(data_sheet) == TimeOfDay::Sunrise);
	}
}



TEST(max_TD,optional)
{
	{
		DataSheet ds;		
		auto m = max_TD(ds);
		ASSERT_FALSE(m.has_value());
	}

	{
		DataSheet ds;
		ds.thermometers_start = {{1,1,1,1}};		
		auto m = max_TD(ds);
		ASSERT_TRUE(*m == 0);
	}

	{
		DataSheet ds;		
		ds.thermometers_end = {{2,0,4,2}};
		auto m = max_TD(ds);
		ASSERT_TRUE(*m == 4);
	}


	{
		DataSheet ds;
		ds.thermometers_start = {{2,0,4,2}};
		ds.thermometers_end = {{1,1,1,1}};
		auto m = max_TD(ds);
		ASSERT_TRUE(*m == 4);
	}
}




TEST(mean_T,optional)
{
	{
		DataSheet ds;
		auto m = mean_T(ds);
		ASSERT_FALSE(m.has_value());
	}

	{
		DataSheet ds;
		ds.thermometers_start = {{2,0,4,2}};		
		auto m = mean_T(ds);
		ASSERT_TRUE(*m == 2);
	}

	{
		DataSheet ds;		
		ds.thermometers_end = {{2,0,4,2}};
		auto m = mean_T(ds);
		ASSERT_TRUE(*m == 2);
	}

	{
		DataSheet ds;
		ds.thermometers_start = {{2,0,4,2}};
		ds.thermometers_end = {{1,1,1,1}};
		auto m = mean_T(ds);
		ASSERT_TRUE(*m == 1.5);
	}

}



TEST(DataSheetStats,test)
{
	{
		DataSheet data_sheet;
		auto stats = data_sheet_stats(data_sheet);
		ASSERT_FALSE(stats.max_TD);
		ASSERT_FALSE(stats.mean_dT);
		ASSERT_FALSE(stats.mean_T);
		//ASSERT_TRUE(stats.abs_drift==0);
		//ASSERT_TRUE(stats.drift==0);
		ASSERT_TRUE(stats.number_of_adjust==0);
		ASSERT_TRUE(stats.number_of_minus==0);
		ASSERT_TRUE(stats.number_of_plus==0);
	}

	{
		DataSheet data_sheet;

		data_sheet.thermometers_start_time = {13,00};
		data_sheet.thermometers_start = {1,3,3,1};
		data_sheet.thermometers_end_time = {13,15};
		data_sheet.thermometers_end = {2,5,5,2};

		{
			DataSheet::Turn turn;
			turn.distances = {0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0 ,0};
			data_sheet.turns.push_back(turn);
		}
		{
			DataSheet::Turn turn;
			turn.bad=true;
			turn.distances = {0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0 ,0};
			data_sheet.turns.push_back(turn);
		}
		{
			DataSheet::Turn turn;
			turn.invert=true;
			turn.distances = {0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0 ,0};
			data_sheet.turns.push_back(turn);
		}
		{
			DataSheet::Turn turn;
			turn.reverse=true;
			turn.distances = {0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0 ,0};
			data_sheet.turns.push_back(turn);
		}
		{ // adjust
			DataSheet::Turn turn;
			turn.adjust = true;
			turn.minus=true;
			turn.distances = {0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0 ,0};
			data_sheet.turns.push_back(turn);
		}
		{
			DataSheet::Turn turn;
			turn.cancel=true;
			turn.distances = {0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0 ,0};
			data_sheet.turns.push_back(turn);
		}
		{
			DataSheet::Turn turn;
			turn.distances = {0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0 ,0};
			data_sheet.turns.push_back(turn);
		}
		{
			DataSheet::Turn turn;
			turn.plus=true;
			turn.distances = {0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0 ,0};
			data_sheet.turns.push_back(turn);
		}

		auto stats = data_sheet_stats(data_sheet);
		ASSERT_TRUE(stats.abs_drift == 0);
		ASSERT_TRUE(stats.drift == 0);
		ASSERT_TRUE(stats.number_of_adjust == 1);
		ASSERT_TRUE(stats.number_of_plus == 2);
		ASSERT_TRUE(stats.number_of_minus == 4);
		ASSERT_TRUE(*stats.max_TD == 3);
		ASSERT_TRUE(*stats.mean_T == 2.75);
		ASSERT_TRUE(*stats.mean_dT == 1.5);
	}

	{
		DataSheet data_sheet;
		data_sheet.thermometers_start = {1,3,3,1};
		auto stats = data_sheet_stats(data_sheet);

		ASSERT_TRUE(*stats.max_TD == 2);
		ASSERT_TRUE(*stats.mean_T == 2);
		ASSERT_FALSE(stats.mean_dT);
	}

	{
		DataSheet data_sheet;
		data_sheet.thermometers_end = {1,3,3,1};
		auto stats = data_sheet_stats(data_sheet);

		ASSERT_TRUE(*stats.max_TD == 2);
		ASSERT_TRUE(*stats.mean_T == 2);
		ASSERT_FALSE(stats.mean_dT);
	}
}




TEST(epoch,unsupported_month)
{
	{
		DataSheet data_sheet;
		data_sheet.date = {1925,1,5};
		try {
			epoch(data_sheet);
			FAIL();
		}
		catch(std::runtime_error&) {			
		}					
	}
	
	{
		DataSheet data_sheet;
		data_sheet.date = {1925,5,5};
		try {
			epoch(data_sheet);
			FAIL();
		}
		catch(std::runtime_error&) {			
		}					
	}
	
	{
		DataSheet data_sheet;
		data_sheet.date = {1925,6,5};
		try {
			epoch(data_sheet);
			FAIL();
		}
		catch(std::runtime_error&) {			
		}					
	}
	
	{
		DataSheet data_sheet;
		data_sheet.date = {1925,10,5};
		try {
			epoch(data_sheet);
			FAIL();
		}
		catch(std::runtime_error&) {			
		}					
	}
	
	{
		DataSheet data_sheet;
		data_sheet.date = {1925,11,5};
		try {
			epoch(data_sheet);
			FAIL();
		}
		catch(std::runtime_error&) {			
		}					
	}
	
	{
		DataSheet data_sheet;
		data_sheet.date = {1925,12,5};
		try {
			epoch(data_sheet);
			FAIL();
		}
		catch(std::runtime_error&) {			
		}					
	}
}


TEST(epoch,test)
{
	{
		DataSheet data_sheet;
		data_sheet.date = {1925,2,5};
		auto e = epoch(data_sheet);
		ASSERT_TRUE(e == Epoch::Feb);
	}
	
	{
		DataSheet data_sheet;
		data_sheet.date = {1925,3,5};
		auto e = epoch(data_sheet);
		ASSERT_TRUE(e == Epoch::Apr);
	}
	
	{
		DataSheet data_sheet;
		data_sheet.date = {1925,4,5};
		auto e = epoch(data_sheet);
		ASSERT_TRUE(e == Epoch::Apr);
	}
	
	{
		DataSheet data_sheet;
		data_sheet.date = {1925,7,5};
		auto e = epoch(data_sheet);
		ASSERT_TRUE(e == Epoch::Aug);
	}
	
	{
		DataSheet data_sheet;
		data_sheet.date = {1925,8,5};
		auto e = epoch(data_sheet);
		ASSERT_TRUE(e == Epoch::Aug);
	}
	
	{
		DataSheet data_sheet;
		data_sheet.date = {1925,9,5};
		auto e = epoch(data_sheet);
		ASSERT_TRUE(e == Epoch::Sep);
	}
}



TEST(load_data_sheet_csv,header)
{
	const std::string name = "dcm1.csv";
	{
		Options options;		
		auto data_sheet = load_data_sheet_csv(name,options);
		
		ASSERT_TRUE(data_sheet.no == 49);
		ASSERT_TRUE(data_sheet.date == DataSheet::Date({1925,9,17}));
		ASSERT_TRUE(data_sheet.mean_observation_time == DataSheet::Time({19,30}));
		ASSERT_TRUE(data_sheet.local_mean_time.has_value());
		ASSERT_TRUE(*data_sheet.local_mean_time == DataSheet::Time({19,38}));
		ASSERT_TRUE(data_sheet.sidereal_time == DataSheet::Time({19,23}));
		ASSERT_TRUE(data_sheet.weight.has_value());
		ASSERT_TRUE(*data_sheet.weight == 8);
		ASSERT_FALSE(data_sheet.fringes.has_value());
		ASSERT_TRUE(data_sheet.sign_correct);
		ASSERT_TRUE(data_sheet.comment == "Test file");
		ASSERT_TRUE(data_sheet.thermometers_start_time.has_value());
		ASSERT_TRUE(*data_sheet.thermometers_start_time == DataSheet::Time({19,18}));
		ASSERT_TRUE(data_sheet.thermometers_start.has_value());
		ASSERT_TRUE(data_sheet.thermometers_start->N == 15.9);
		ASSERT_TRUE(data_sheet.thermometers_start->E == 15.8);
		ASSERT_TRUE(data_sheet.thermometers_start->S == 15.7);
		ASSERT_TRUE(data_sheet.thermometers_start->W == 15.6);
		ASSERT_TRUE(data_sheet.weather_clouds);
		ASSERT_FALSE(data_sheet.weather_rain);
		ASSERT_FALSE(data_sheet.weather_wind);
		ASSERT_TRUE(data_sheet.weather_obscuration == DataSheet::WeatherObscurations::Clear);
		ASSERT_TRUE(data_sheet.thermometers_end_time.has_value());
		ASSERT_TRUE(*data_sheet.thermometers_end_time == DataSheet::Time({19,41}));
		ASSERT_TRUE(data_sheet.thermometers_end->N == 15.98);
		ASSERT_TRUE(data_sheet.thermometers_end->E == 15.68);
		ASSERT_TRUE(data_sheet.thermometers_end->S == 15.8);
		ASSERT_TRUE(data_sheet.thermometers_end->W == 15.9);
		ASSERT_TRUE(data_sheet.start_time == DataSheet::Time({19,20}));
		ASSERT_TRUE(data_sheet.end_time == DataSheet::Time({19,41}));
		ASSERT_FALSE(data_sheet.reverse);
		ASSERT_FALSE(data_sheet.inverted);
		ASSERT_FALSE(data_sheet.desk_in_sw);
		ASSERT_FALSE(data_sheet.corrugated_paper);
	}
}



TEST(load_data_sheet_csv,codes)
{	
	const std::string name = "dcm1.csv";
	
	{
		Options options;		
		auto data_sheet = load_data_sheet_csv(name,options);
		
		ASSERT_TRUE(data_sheet.turns.size() == 20u);
		ASSERT_TRUE(data_sheet.turns.at(0).reverse);
		ASSERT_TRUE(data_sheet.turns.at(1).cancel);
		ASSERT_TRUE(!data_sheet.turns.at(2).reverse);
		ASSERT_TRUE(!data_sheet.turns.at(3).cancel);
		ASSERT_TRUE(data_sheet.turns.at(4).bad);
		ASSERT_TRUE(data_sheet.turns.at(5).invert);
		ASSERT_TRUE(data_sheet.turns.at(6).bad && data_sheet.turns.at(6).cancel);
		ASSERT_TRUE(data_sheet.turns.at(7).invert && data_sheet.turns.at(7).reverse);
		ASSERT_TRUE(data_sheet.turns.at(8).bad && !data_sheet.turns.at(8).cancel);
		ASSERT_TRUE(data_sheet.turns.at(9).invert && !data_sheet.turns.at(9).reverse);
		
	}
	
	{
		Options options;
		options.ignore_bad = true;
		options.ignore_invert = true;		
		options.ignore_cancel_disabled = true;
		options.ignore_reverse_disabled = true;		
		auto data_sheet = load_data_sheet_csv(name,options);
		
		ASSERT_TRUE(data_sheet.turns.size() == 20u);
		ASSERT_TRUE(data_sheet.turns.at(0).reverse);
		ASSERT_TRUE(data_sheet.turns.at(1).cancel);
		ASSERT_TRUE(data_sheet.turns.at(2).reverse);
		ASSERT_TRUE(data_sheet.turns.at(3).cancel);
		ASSERT_TRUE(!data_sheet.turns.at(4).bad);
		ASSERT_TRUE(!data_sheet.turns.at(5).invert);
		ASSERT_TRUE(!data_sheet.turns.at(6).bad && data_sheet.turns.at(6).cancel);
		ASSERT_TRUE(!data_sheet.turns.at(7).invert && data_sheet.turns.at(7).reverse);
		ASSERT_TRUE(!data_sheet.turns.at(8).bad && data_sheet.turns.at(8).cancel);
		ASSERT_TRUE(!data_sheet.turns.at(9).invert && data_sheet.turns.at(9).reverse);
		
	}
	
	
	{
		Options options;
		options.ignore_bad = true;
		options.ignore_cancel = true;
		options.ignore_invert = true;
		options.ignore_reverse = true;
		options.ignore_cancel_disabled = true;
		options.ignore_reverse_disabled = true;
		auto data_sheet = load_data_sheet_csv(name,options);
		
		ASSERT_TRUE(data_sheet.turns.size() == 20u);
		ASSERT_TRUE(!data_sheet.turns.at(0).reverse);
		ASSERT_TRUE(!data_sheet.turns.at(1).cancel);
		ASSERT_TRUE(!data_sheet.turns.at(2).reverse);
		ASSERT_TRUE(!data_sheet.turns.at(3).cancel);
		ASSERT_TRUE(!data_sheet.turns.at(4).bad);
		ASSERT_TRUE(!data_sheet.turns.at(5).invert);
		ASSERT_TRUE(!data_sheet.turns.at(6).bad && !data_sheet.turns.at(6).cancel);
		ASSERT_TRUE(!data_sheet.turns.at(7).invert && !data_sheet.turns.at(7).reverse);
		ASSERT_TRUE(!data_sheet.turns.at(8).bad && !data_sheet.turns.at(8).cancel);
		ASSERT_TRUE(!data_sheet.turns.at(9).invert && !data_sheet.turns.at(9).reverse);
		
	}
}



int main(int ,char** )
{
	RUN_ALL_TESTS();

	return EXIT_SUCCESS;
}
