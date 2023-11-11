cmake_minimum_required(VERSION 2.8.3)
project(inekf)

# no warnings 
add_definitions(-w)

option(USE_CPP11 "Use flag -std=c++11" ON)
option(USE_MUTEX "Use mutex in code" OFF)

message("USE_CPP11: " ${USE_CPP11})

if (USE_CPP11)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
endif ()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -DEIGEN_NO_DEBUG -march=native -Wl,--no-as-needed")
SET(CMAKE_CXX_COMPILER /usr/bin/g++)

message("CMAKE_CXX_FLAGS: " ${CMAKE_CXX_FLAGS})
unset(USE_CPP11)
unset(USE_MUTEX)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)

include_directories("${PROJECT_SOURCE_DIR}/include/")

# Threading
find_package(Threads)

# Boost
find_package(Boost 1.57 REQUIRED COMPONENTS system)
include_directories(${Boost_INCLUDE_DIR})
message("Boost_INCLUDE_DIR: " ${Boost_INCLUDE_DIR})

# Eigen 
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})
message("EIGEN3_INCLUDE_DIR: " ${EIGEN3_INCLUDE_DIR})

# OpenMP
#find_package(OpenMP REQUIRED)
#if (OPENMP_FOUND)
#	message("OpenMP: Found")
#	# cmake setting
#	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
#    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${OpenMP_SHARED_LINKER_FLAGS}")
#endif ()

file(GLOB src_files
		"${PROJECT_SOURCE_DIR}/src/*.cpp"
		)

add_executable(landmarks ${PROJECT_SOURCE_DIR}/src/examples/landmarks.cpp ${src_files})
add_executable(kinematics ${PROJECT_SOURCE_DIR}/src/examples/kinematics.cpp ${src_files})
add_executable(propagation_speed ${PROJECT_SOURCE_DIR}/src/tests/propagation_speed.cpp ${src_files})
add_executable(correction_speed ${PROJECT_SOURCE_DIR}/src/tests/correction_speed.cpp ${src_files})
add_executable(large_state_speed ${PROJECT_SOURCE_DIR}/src/tests/large_state_speed.cpp ${src_files})
add_executable(left_vs_right_error_dynamics ${PROJECT_SOURCE_DIR}/src/tests/left_vs_right_error_dynamics.cpp ${src_files})

target_link_libraries(landmarks ${Boost_LIBRARIES})
target_link_libraries(kinematics ${Boost_LIBRARIES})
target_link_libraries(propagation_speed ${Boost_LIBRARIES})
target_link_libraries(correction_speed ${Boost_LIBRARIES})
target_link_libraries(large_state_speed ${Boost_LIBRARIES})
target_link_libraries(left_vs_right_error_dynamics ${Boost_LIBRARIES})

add_library(inekf SHARED
		"${PROJECT_SOURCE_DIR}/src/inekf.cpp"
		"${PROJECT_SOURCE_DIR}/src/state.cpp"
		"${PROJECT_SOURCE_DIR}/src/matrix_lie_group.cpp"
		"${PROJECT_SOURCE_DIR}/src/noise.cpp"
		"${PROJECT_SOURCE_DIR}/src/observation.cpp"
		)

set_target_properties(inekf PROPERTIES LIBRARY_OUTPUT_NAME inekf)

# external library
export(PACKAGE inekf)

#set_target_properties(inekf PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
get_property(inekf_include_dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
get_property(inekf_library_dirs TARGET inekf PROPERTY LIBRARY_OUTPUT_DIRECTORY)
get_property(inekf_libraries TARGET inekf PROPERTY LIBRARY_OUTPUT_NAME)

message("inekf_include_dirs: " ${inekf_include_dirs})
message("inekf_library_dirs: " ${inekf_library_dirs})
message("inekf_libraries: " ${inekf_libraries})

configure_file(inekfConfig.cmake.in
		"${PROJECT_BINARY_DIR}/inekfConfig.cmake" @ONLY
		)

message("PROJECT_BINARY_DIR: " ${PROJECT_BINARY_DIR})
