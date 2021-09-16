

#include "test.h"

#include "../Filter.h"
#include "../DataSheet.h"

using namespace aether;


DataSheet create_test_sheet()
{
	DataSheet data_sheet;
	data_sheet.date = {1925,2,1};
	data_sheet.no = 14;
	data_sheet.start_time = {12,30};
	data_sheet.end_time = {12,46};
	data_sheet.mean_observation_time = {12,38};
	data_sheet.local_mean_time = {12,46};
	data_sheet.sidereal_time = {6,38}; // amplitude 0.04
	data_sheet.fringes = 6;
	data_sheet.weight = 8;
	data_sheet.thermometers_start = {1,2,3,2};
	data_sheet.thermometers_end = {1,2,3.5,2};
	data_sheet.desk_in_sw = false;


	{
		DataSheet::Turn turn;
		turn.distances = {0,1,2,3,2,1,0,-1, -2,-3,-2,-1,0,1,2,3, 2};
		data_sheet.turns.push_back(turn);
	}

	return data_sheet;
}


TEST(Filter,test)
{
	{
		Filter::Interval interval;
		ASSERT_TRUE(std::isinf(interval.min));
		ASSERT_TRUE(std::isinf(interval.max));
	}
	{
		DataSheet data_sheet;
		Options options;
		Filter filter;
		ASSERT_TRUE(selected(data_sheet,options,filter));
	}

	{
		DataSheet data_sheet = create_test_sheet();
		Options options;
		Filter filter;
		ASSERT_TRUE(selected(data_sheet,options,filter));
	}



}




TEST(Filter,year)
{
	{
		DataSheet data_sheet = create_test_sheet();
		Options options;
		Filter filter;
		filter.year.push_back(Filter::Interval{1925,INFINITY});
		ASSERT_TRUE(selected(data_sheet,options,filter));
	}

	{
		DataSheet data_sheet = create_test_sheet();
		Options options;
		Filter filter;
		filter.year.push_back(Filter::Interval{-INFINITY,1924});
		ASSERT_FALSE(selected(data_sheet,options,filter));
	}
}



TEST(Filter,month)
{
	{
		DataSheet data_sheet = create_test_sheet();
		Options options;
		Filter filter;
		filter.month.push_back(Filter::Interval{1,3});
		ASSERT_TRUE(selected(data_sheet,options,filter));
	}

	{
		DataSheet data_sheet = create_test_sheet();
		Options options;
		Filter filter;
		filter.month.push_back(Filter::Interval{3,4});
		ASSERT_FALSE(selected(data_sheet,options,filter));
	}
}



TEST(Filter,day)
{
	{
		DataSheet data_sheet = create_test_sheet();
		Options options;
		Filter filter;
		filter.day.push_back(Filter::Interval{1,3});
		ASSERT_TRUE(selected(data_sheet,options,filter));
	}

	{
		DataSheet data_sheet = create_test_sheet();
		Options options;
		Filter filter;
		filter.day.push_back(Filter::Interval{3,4});
		ASSERT_FALSE(selected(data_sheet,options,filter));
	}
}




TEST(Filter,no)
{
	{
		DataSheet data_sheet = create_test_sheet();
		Options options;
		Filter filter;
		filter.no.push_back(Filter::Interval{13,15});
		ASSERT_TRUE(selected(data_sheet,options,filter));
	}

	{
		DataSheet data_sheet = create_test_sheet();
		Options options;
		Filter filter;
		filter.no.push_back(Filter::Interval{23,24});
		ASSERT_FALSE(selected(data_sheet,options,filter));
	}
}



TEST(Filter,time)
{
	{
		DataSheet data_sheet = create_test_sheet();
		Options options;
		Filter filter;
		filter.time.push_back(Filter::Interval{12.5,13});
		ASSERT_TRUE(selected(data_sheet,options,filter));
	}

	{
		DataSheet data_sheet = create_test_sheet();
		Options options;
		Filter filter;
		filter.time.push_back(Filter::Interval{13,14});
		ASSERT_FALSE(selected(data_sheet,options,filter));
	}
}



TEST(Filter,sidereal_time)
{
	{
		DataSheet data_sheet = create_test_sheet();
		Options options;
		Filter filter;
		filter.sidereal_time.push_back(Filter::Interval{6,7});
		ASSERT_TRUE(selected(data_sheet,options,filter));
	}

	{
		DataSheet data_sheet = create_test_sheet();
		Options options;
		Filter filter;
		filter.sidereal_time.push_back(Filter::Interval{7,8});
		ASSERT_FALSE(selected(data_sheet,options,filter));
	}
}



TEST(Filter,adjust)
{
	{
		DataSheet data_sheet = create_test_sheet();
		Options options;
		Filter filter;
		filter.adjust.push_back(Filter::Interval{0,3});
		ASSERT_TRUE(selected(data_sheet,options,filter));
	}

	{
		DataSheet data_sheet = create_test_sheet();
		Options options;
		Filter filter;
		filter.adjust.push_back(Filter::Interval{1,3});
		ASSERT_FALSE(selected(data_sheet,options,filter));
	}
}




TEST(Filter,fringes)
{
	{
		DataSheet data_sheet = create_test_sheet();
		Options options;
		Filter filter;
		filter.fringes.push_back(Filter::Interval{6,7});
		ASSERT_TRUE(selected(data_sheet,options,filter));
	}

	{
		DataSheet data_sheet = create_test_sheet();
		Options options;
		Filter filter;
		filter.fringes.push_back(Filter::Interval{1,3});
		ASSERT_FALSE(selected(data_sheet,options,filter));
	}
}




TEST(Filter,weight)
{
	{
		DataSheet data_sheet = create_test_sheet();
		Options options;
		Filter filter;
		filter.weight.push_back(Filter::Interval{7,8});
		ASSERT_TRUE(selected(data_sheet,options,filter));
	}

	{
		DataSheet data_sheet = create_test_sheet();
		Options options;
		Filter filter;
		filter.weight.push_back(Filter::Interval{8.1,9});
		ASSERT_FALSE(selected(data_sheet,options,filter));
	}
}



TEST(Filter,T)
{
	{
		DataSheet data_sheet = create_test_sheet();
		Options options;
		Filter filter;
		filter.T.push_back(Filter::Interval{1,4});
		ASSERT_TRUE(selected(data_sheet,options,filter));
	}

	{
		DataSheet data_sheet = create_test_sheet();
		Options options;
		Filter filter;
		filter.T.push_back(Filter::Interval{1,3});
		ASSERT_FALSE(selected(data_sheet,options,filter));
	}
}




TEST(Filter,TD)
{
	{
		DataSheet data_sheet = create_test_sheet();
		Options options;
		Filter filter;
		filter.TD.push_back(Filter::Interval{2,3});
		ASSERT_TRUE(selected(data_sheet,options,filter));
	}

	{
		DataSheet data_sheet = create_test_sheet();
		Options options;
		Filter filter;
		filter.TD.push_back(Filter::Interval{3,9});
		ASSERT_FALSE(selected(data_sheet,options,filter));
	}
}




TEST(Filter,dT)
{
	{
		DataSheet data_sheet = create_test_sheet();
		Options options;
		Filter filter;
		filter.dT.push_back(Filter::Interval{0,1});
		ASSERT_TRUE(selected(data_sheet,options,filter));
	}

	{
		DataSheet data_sheet = create_test_sheet();
		Options options;
		Filter filter;
		filter.dT.push_back(Filter::Interval{0.6,1});
		ASSERT_FALSE(selected(data_sheet,options,filter));
	}
}




TEST(Filter,mean_dT)
{
	{
		DataSheet data_sheet = create_test_sheet();
		Options options;
		Filter filter;
		filter.mean_dT.push_back(Filter::Interval{0.1,0.2});
		ASSERT_TRUE(selected(data_sheet,options,filter));
	}

	{
		DataSheet data_sheet = create_test_sheet();
		Options options;
		Filter filter;
		filter.mean_dT.push_back(Filter::Interval{0,0.1});
		ASSERT_FALSE(selected(data_sheet,options,filter));
	}
}


TEST(Filter,amplitude)
{
	{
		DataSheet data_sheet = create_test_sheet(); // amplitude = 0.3
		Options options;
		Filter filter;
		filter.amplitude.push_back(Filter::Interval{0.1,0.2});
		ASSERT_TRUE(selected(data_sheet,options,filter));
	}

	{
		DataSheet data_sheet = create_test_sheet();
		Options options;
		Filter filter;
		filter.amplitude.push_back(Filter::Interval{0.1,0.149});
		ASSERT_FALSE(selected(data_sheet,options,filter));
	}
}



TEST(Filter,drift)
{
	{
		DataSheet data_sheet = create_test_sheet();
		Options options;
		Filter filter;
		filter.drift.push_back(Filter::Interval{1,3});
		ASSERT_TRUE(selected(data_sheet,options,filter));
	}

	{
		DataSheet data_sheet = create_test_sheet();
		Options options;
		Filter filter;
		filter.drift.push_back(Filter::Interval{3,4});
		ASSERT_FALSE(selected(data_sheet,options,filter));
	}
}




TEST(Filter,abs_drift)
{
	{
		DataSheet data_sheet = create_test_sheet();
		Options options;
		Filter filter;
		filter.abs_drift.push_back(Filter::Interval{1,3});
		ASSERT_TRUE(selected(data_sheet,options,filter));
	}

	{
		DataSheet data_sheet = create_test_sheet();
		Options options;
		Filter filter;
		filter.abs_drift.push_back(Filter::Interval{3,4});
		ASSERT_FALSE(selected(data_sheet,options,filter));
	}
}







TEST(Filter,uncertainty)
{
	DataSheet data_sheet;
	for (int i=0;i<2;i++) {
		DataSheet::Turn turn;
		turn.distances = {0,1,2,3,2,1,0,-1, -2,-3,-2,-1,0,1,2,3, 2};
		data_sheet.turns.push_back(turn);
	}

	{
		Options options;
		Filter filter;
		filter.uncertainty.push_back(Filter::Interval{0,1});
		ASSERT_TRUE(selected(data_sheet,options,filter));
	}

	{
		Options options;
		Filter filter;
		filter.uncertainty.push_back(Filter::Interval{1,2});
		ASSERT_FALSE(selected(data_sheet,options,filter));
	}
}




TEST(Filter,theory_amplitude)
{
	{
		DataSheet data_sheet = create_test_sheet();
		Options options;
		Filter filter;
		filter.theory_amplitude.push_back(Filter::Interval{0.017,0.023});
		ASSERT_TRUE(selected(data_sheet,options,filter));
	}

	{
		DataSheet data_sheet = create_test_sheet();
		Options options;
		Filter filter;
		filter.theory_amplitude.push_back(Filter::Interval{0,0.005});
		ASSERT_FALSE(selected(data_sheet,options,filter));
	}
}




TEST(Filter,nw)
{
	{
		DataSheet data_sheet;
		data_sheet.desk_in_sw = false;
		Options options;
		Filter filter;
		filter.nw = true;
		ASSERT_TRUE(selected(data_sheet,options,filter));
	}

	{
		DataSheet data_sheet;
		data_sheet.desk_in_sw = true;
		Options options;
		Filter filter;
		filter.nw = true;
		ASSERT_FALSE(selected(data_sheet,options,filter));
	}
}



TEST(Filter,sw)
{
	{
		DataSheet data_sheet;
		data_sheet.desk_in_sw = true;
		Options options;
		Filter filter;
		filter.sw = true;
		ASSERT_TRUE(selected(data_sheet,options,filter));
	}

	{
		DataSheet data_sheet;
		data_sheet.desk_in_sw = false;
		Options options;
		Filter filter;
		filter.sw = true;
		ASSERT_FALSE(selected(data_sheet,options,filter));
	}
}



int main(int ,char** )
{
	RUN_ALL_TESTS();

	return EXIT_SUCCESS;
}
