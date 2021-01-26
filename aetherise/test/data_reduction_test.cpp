
#include "test.h"

#include "../data_reduction.h"

using namespace aether;


bool approximate(const std::array<double,17>& a,const std::array<double,17>& b,double g)
{
	for (int i=0;i<17;i++) {
		if (!test::approximate(a.at(i),b.at(i),g))
			return false;
	}
	return true;
}

const std::array<float,17> test_distances = {0,1,2,3,4,3,2,1, 0,-1,-2,-3,-4,-3,-2,-1, 0};


DataSheet create_data_sheet()
{

	DataSheet data_sheet;

	DataSheet::Turn turn;
	turn.distances = test_distances;
	data_sheet.turns.push_back(turn);

	//turn.distances = test_distances;
	int drift=0;
	for (auto& d : turn.distances)
		d = d+3 + drift++;
	data_sheet.turns.push_back(turn);

	turn.distances = test_distances;
	drift=0;
	for (auto& d : turn.distances)
		d = d+8 + drift--;
	data_sheet.turns.push_back(turn);

	return data_sheet;
}

DataSheet create_data_sheet_with_codes()
{

	DataSheet data_sheet;

	for (int i=0;i<10;i++) {
		DataSheet::Turn turn;
		turn.distances = test_distances;
		turn.distances.at(0)=i;
		data_sheet.turns.push_back(turn);
	}
	data_sheet.turns.at(0).bad = true;
	data_sheet.turns.at(1).invert = true;
	data_sheet.turns.at(2).cancel = true;
	data_sheet.turns.at(3).reverse = true;
	
	return data_sheet;
}



TEST(MillersReduction,test)
{
	{
		DataSheet data_sheet = create_data_sheet();
		Options options;		

		ReducedData reduced_data = reduce_data(data_sheet,options);
		std::array<double,17> expected = {0,0.1,0.2,0.3,0.4,0.3,0.2,0.1,0, -0.1,-0.2,-0.3,-0.4,-0.3,-0.2,-0.1, 0};
		for (auto& e : expected)
			e*=0.5;
		ASSERT_TRUE(approximate(reduced_data.displacements,expected,0.001));
		std::array<double,17> expected_u {};
		ASSERT_TRUE(reduced_data.uncertainties==expected_u);
	}
}






TEST(reduce_to_single_period,test)
{
	{
		ReducedData reduced_data;
		reduced_data.displacements = {1,1,1,1,1,1,1,1, 3,3,3,3,3,3,3,3, 0};
		reduced_data.uncertainties = {4,4,4,4,4,4,4,4, 3,3,3,3,3,3,3,3, 4};

		std::array<double,17> expected_d {2,2,2,2,2,2,2,2, 1.5,2,2,2,2,2,2,2, 1.5};
		reduce_to_single_period(reduced_data);

		ASSERT_TRUE(reduced_data.displacements==expected_d);
		for (int i=0;i<17;i++) {
			ASSERT_APPROX(reduced_data.uncertainties.at(i),(5.0/2.0),0.001);
		}

	}

	{
		ReducedData reduced_data;
		reduced_data.displacements = {0, 1, 2, 3, 4, 3, 2, 1,
									  0,-1,-2,-3,-4,-3,-2,-1, 0};
		reduced_data.uncertainties = {1,1,1,1,1,1,1,1, 1,1,1,1,1,1,1,1,1};

		std::array<double,17> expected_d {0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0};
		reduce_to_single_period(reduced_data);

		ASSERT_TRUE(reduced_data.displacements==expected_d);
		for (int i=0;i<17;i++) {
			ASSERT_APPROX(reduced_data.uncertainties.at(i),std::sqrt(2)/2,0.001);
		}

	}

	{
		ReducedData reduced_data;
		reduced_data.displacements = {-1, 0, 1, 2, 3, 4, 3, 2,
									   1, 0,-1,-2,-3,-4,-3,-2,-1};
		reduced_data.uncertainties = {1,1,1,1,1,1,1,1, 1,1,1,1,1,1,1,1,1};

		std::array<double,17> expected_d {0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0};
		reduce_to_single_period(reduced_data);

		ASSERT_TRUE(reduced_data.displacements==expected_d);
		for (int i=0;i<17;i++) {
			ASSERT_APPROX(reduced_data.uncertainties.at(i),std::sqrt(2)/2,0.001);
		}

	}

	{
		ReducedData reduced_data;
		reduced_data.displacements = { 0, 1, 2, 1, 0, -1, -2, -1,
									   0, 1, 2, 1, 0, -1, -2, -1, 0};
		reduced_data.uncertainties = {1,1,1,1,1,1,1,1, 1,1,1,1,1,1,1,1,1};

		std::array<double,17> expected_d {0,1,2,1,0,-1,-2,-1, 0,1,2,1,0,-1,-2,-1,0};
		reduce_to_single_period(reduced_data);

		ASSERT_TRUE(reduced_data.displacements==expected_d);
		for (int i=0;i<17;i++) {
			ASSERT_APPROX(reduced_data.uncertainties.at(i),std::sqrt(2)/2,0.001);
		}

	}
}



TEST(selected_and_tranformed_turns,test)
{
	{
		Options options;
		DataSheet data_sheet = create_data_sheet_with_codes();
			
		
		std::vector<std::array<float,17>> sat;
		selected_and_transformed_turns(data_sheet,options,[&sat](int ,const DataSheet::Turn& ,const std::array<float,17>& d){
			sat.push_back(d);
		});
		
		ASSERT_TRUE(sat.size()==8u);		
		ASSERT_TRUE(sat.at(0).at(0)==-1);		
		ASSERT_TRUE(sat.at(1).at(0)==-3);
		ASSERT_TRUE(sat.at(2).at(0)== 4);
	}
	
	{
		Options options;		
		DataSheet data_sheet = create_data_sheet_with_codes();
		data_sheet.reverse = true;
			
		
		std::vector<std::array<float,17>> sat;
		selected_and_transformed_turns(data_sheet,options,[&sat](int ,const DataSheet::Turn& ,const std::array<float,17>& d){
			sat.push_back(d);
		});
		
		ASSERT_TRUE(sat.size()==8u);		
		ASSERT_TRUE(sat.at(0).at(0)==1);		
		ASSERT_TRUE(sat.at(1).at(0)==3);
		ASSERT_TRUE(sat.at(2).at(0)==-4);
	}
	
	{
		Options options;		
		options.invert_data = true;
		DataSheet data_sheet = create_data_sheet_with_codes();		
			
		
		std::vector<std::array<float,17>> sat;
		selected_and_transformed_turns(data_sheet,options,[&sat](int ,const DataSheet::Turn& ,const std::array<float,17>& d){
			sat.push_back(d);
		});
		
		ASSERT_TRUE(sat.size()==8u);		
		ASSERT_TRUE(sat.at(0).at(0)==1);		
		ASSERT_TRUE(sat.at(1).at(0)==3);
		ASSERT_TRUE(sat.at(2).at(0)==-4);
	}
	
	{
		Options options;
		options.invert_data = true;
		DataSheet data_sheet = create_data_sheet_with_codes();
		data_sheet.reverse = true;
			
		
		std::vector<std::array<float,17>> sat;
		selected_and_transformed_turns(data_sheet,options,[&sat](int ,const DataSheet::Turn& ,const std::array<float,17>& d){
			sat.push_back(d);
		});
		
		ASSERT_TRUE(sat.size()==8u);		
		ASSERT_TRUE(sat.at(0).at(0)==-1);		
		ASSERT_TRUE(sat.at(1).at(0)==-3);
		ASSERT_TRUE(sat.at(2).at(0)==4);
	}
}



int main(int ,char** )
{
	RUN_ALL_TESTS();

	return EXIT_SUCCESS;
}
