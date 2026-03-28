# ########################################################################################
# Use the code in this section for standalone projects
cmake_minimum_required(VERSION 3.12 FATAL_ERROR)
project(RPC)
if(NOT TARGET Garfield::Garfield)
  find_package(Garfield REQUIRED)
endif()
# ########################################################################################

# ---Build executable------------------------------------------------------------
add_executable(Weight Weight.C)
add_executable(Weight_V2 Weight_V2.C)
add_executable(Test Test.C)
add_executable(Test_10kV Test_10kV.C)
add_executable(Test_10kV_500Step Test_10kV_500Step.C)
target_link_libraries(Test_10kV_500Step Garfield::Garfield)
target_link_libraries(Weight Garfield::Garfield)
target_link_libraries(Weight_V2 Garfield::Garfield)
target_link_libraries(Test_10kV Garfield::Garfield)
target_link_libraries(Test Garfield::Garfield)
# ---Copy all files locally to the build directory-------------------------------
foreach(_file c2h2f4_ic4h10_sf6.gas TimingRPCGas.gas)
  configure_file(${_file} ${_file} COPYONLY)
endforeach()
