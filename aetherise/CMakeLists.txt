cmake_minimum_required(VERSION 3.7)

project(aetherise)
set(CMAKE_CXX_STANDARD 11)

option(minuit2 "Enable usage of ROOT::Minuit2 stand alone library")
set(minuit2_include_path "/usr/local/include/Minuit2" CACHE PATH "Path to ROOT::Minuit2 stand alone headers")
set(minuit2_math_lib_file "/usr/local/lib/libMinuit2Math.a" CACHE FILEPATH "Location of libMinuit2Math.a")
set(minuit2_lib_file "/usr/local/lib/libMinuit2.a" CACHE FILEPATH "Location of libMinuit2.a")

if(MSVC)
  add_compile_options(/W4 /WX)
elseif(CMAKE_COMPILER_IS_GNUCXX)
  add_compile_options(-Wall -Wextra -pedantic)
endif()

if(WIN32)
	add_definitions(-DAETHER_WINDOWS)
endif()


message(STATUS "Building with ROOT::Minuit2 stand alone library: " ${minuit2})

if (minuit2)
	# Minuit2 Stand Alone installation

	add_definitions(-DAETHER_MINUIT)

	include_directories(${minuit2_include_path})

	add_library(Minuit2Math STATIC IMPORTED)
	add_library(Minuit2 STATIC IMPORTED)

	set_target_properties(Minuit2Math PROPERTIES IMPORTED_LOCATION ${minuit2_math_lib_file})
	set_target_properties(Minuit2 PROPERTIES IMPORTED_LOCATION ${minuit2_lib_file})
endif()


add_library(aetherise_lib STATIC
	"stdx.h" "stdx.cpp"
	"utils.h" "utils.cpp"
	"mathematics.h" "mathematics.cpp"
	"astro.h" "astro.cpp"
	"DataSheet.h" "DataSheet.cpp"
	"data_reduction.h" "data_reduction.cpp"	
	"Theory.h" "Theory.cpp"
	"physics.h" "physics.cpp"
	"cmd_line.h" "cmd_line.cpp"
	"models.h" "models.cpp"
	"generators.h" "generators.cpp"	
	"mini.h" "mini.cpp"
	"aetherise.h" "aetherise.cpp"
	"Filter.h" "Filter.cpp")

if (minuit2)
	target_link_libraries(aetherise_lib Minuit2 Minuit2Math)
endif()

add_executable(${PROJECT_NAME} "main.cpp")
target_link_libraries(${PROJECT_NAME} aetherise_lib)



enable_testing()
add_executable(calculations_test "test/test.cpp" "test/test.h" "test/calculations.cpp")
target_link_libraries(calculations_test aetherise_lib)
add_test(test_calculations calculations_test)

add_executable(Theory_test "test/test.cpp" "test/test.h" "test/Theory_test.cpp")
target_link_libraries(Theory_test aetherise_lib)
add_test(test_Theory Theory_test)

add_executable(utils_test "test/test.cpp" "test/test.h" "test/utils_test.cpp")
target_link_libraries(utils_test aetherise_lib)
add_test(test_utils utils_test)

add_executable(data_reduction_test "test/test.cpp" "test/test.h" "test/data_reduction_test.cpp")
target_link_libraries(data_reduction_test aetherise_lib)
add_test(test_data_reduction data_reduction_test)

add_executable(DataSheet_test "test/test.cpp" "test/test.h" "test/DataSheet_test.cpp")
target_link_libraries(DataSheet_test aetherise_lib)
add_test(test_DataSheet DataSheet_test)
configure_file("test/dcm1.csv" "dcm1.csv")
set_tests_properties(test_DataSheet PROPERTIES REQUIRED_FILES "dcm1.csv")

add_executable(mathematics_test "test/test.cpp" "test/test.h" "test/mathematics_test.cpp")
target_link_libraries(mathematics_test aetherise_lib)
add_test(test_mathematics mathematics_test)

add_executable(Filter_test "test/test.cpp" "test/test.h" "test/Filter_test.cpp")
target_link_libraries(Filter_test aetherise_lib)
add_test(test_Filter Filter_test)

add_executable(astro_test "test/test.cpp" "test/test.h" "test/astro_test.cpp")
target_link_libraries(astro_test aetherise_lib)
add_test(test_astro astro_test)

add_executable(aetherise_test "test/test.cpp" "test/test.h" "test/aetherise_test.cpp")
target_link_libraries(aetherise_test aetherise_lib)
add_test(test_aetherise aetherise_test) #  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}

add_executable(mini_test "test/test.cpp" "test/test.h" "test/mini_test.cpp")
target_link_libraries(mini_test aetherise_lib)
add_test(test_mini mini_test)

add_executable(physics_test "test/test.cpp" "test/test.h" "test/physics_test.cpp")
target_link_libraries(physics_test aetherise_lib)
add_test(test_physics physics_test)

add_executable(stdx_test "test/test.cpp" "test/test.h" "test/stdx_test.cpp")
target_link_libraries(stdx_test aetherise_lib)
add_test(test_stdx stdx_test)

