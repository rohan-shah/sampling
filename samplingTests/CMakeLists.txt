cmake_minimum_required(VERSION 3.1)

add_executable(samplingTests main.cpp conditionalPoissonSimple.cpp conditionalPoissonDeterministic.cpp conditionalPoissonZeroWeights.cpp conditionalPoissonZeroWeightsAndDeterministic.cpp paretoDeterministic.cpp paretoSimple.cpp paretoZeroWeights.cpp conditionalPoissonDrafting.cpp conditionalPoissonSequential.cpp)

target_include_directories(samplingTests PUBLIC ${Boost_INCLUDE_DIRS} ${PROJECT_SOURCE_DIR})
target_link_libraries(samplingTests PUBLIC ${Boost_LIBRARIES} sampling)
target_compile_definitions(samplingTests PUBLIC BOOST_ALL_NO_LIB BOOST_ALL_DYN_LINK)

if(WIN32)
	set_target_properties(samplingTests PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR} PDB_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
endif()