cmake_minimum_required(VERSION 2.8)

if (NOT DEFINED LEMONADE_INCLUDE_DIR)
message("LEMONADE_INCLUDE_DIR is not provided. If build fails, use -DLEMONADE_INCLUDE_DIR=/path/to/LeMonADE/headers/ or install to default location")
endif()

if (NOT DEFINED LEMONADE_LIBRARY_DIR)
message("LEMONADE_LIBRARY_DIR is not provided. If build fails, use -DLEMONADE_LIBRARY_DIR=/path/to/LeMonADE/lib/ or install to default location")
endif()

include_directories (${LEMONADE_INCLUDE_DIR})
link_directories (${LEMONADE_LIBRARY_DIR})

add_executable(LinearChain LinearChain.cpp)

target_link_libraries(LinearChain LeMonADE )

add_executable(CalcRg2 CalcRg2.cpp)

target_link_libraries(CalcRg2 LeMonADE )

add_executable(CalcRg2Cos CalcRg2Cos.cpp)

target_link_libraries(CalcRg2Cos LeMonADE )

add_executable(CoSolvency CoSolvency.cpp)

target_link_libraries(CoSolvency LeMonADE )

