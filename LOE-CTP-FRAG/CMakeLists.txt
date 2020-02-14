# CMake version
cmake_minimum_required(VERSION 3.12)

# Project declaration
project(LOE-CTP-FRAG)

# Set executable directory
set(EXECUTABLE_OUTPUT_PATH bin)


add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/QC_Tools ${CMAKE_CURRENT_BINARY_DIR}/QC_Tools/QC_Tools_Build)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/Approx_MVBB/ApproxMVBB ${CMAKE_CURRENT_BINARY_DIR}/Approx_MVBB/ApproxMVBB/ApproxMVBB_Build)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

set(BUILD_FLAGS "-Wall -g -std=c++14 -fopenmp -O2 -ip -xCORE-AVX512 -Wl,-rpath,Approx_MVBB/ApproxMVBB/ApproxMVBB_Build/lib Approx_MVBB/ApproxMVBB/ApproxMVBB_Build/lib/libApproxMVBB.so.2.1.0")

#set(BUILD_FLAGS "-Wall -g -std=c++14 -fopenmp -Wl,-rpath,Approx_MVBB/ApproxMVBB/ApproxMVBB_Build/lib Approx_MVBB/ApproxMVBB/ApproxMVBB_Build/lib/libApproxMVBB.so.2.1.0")


# Inclusion of external libraries
include_directories(
	Approx_MVBB/ApproxMVBB/ApproxMVBB_Build/include
	Approx_MVBB/ApproxMVBB/include
	Approx_MVBB/ApproxMVBB/external/GeometryPredicates/include
	Approx_MVBB/ApproxMVBB/external/Diameter/include
	Approx_MVBB/ApproxMVBB/eigen3
	${CMAKE_CURRENT_SOURCE_DIR}
	)
link_directories(
	QC_Tools 
	QC_Tools/QCTools_Build/src/IO 
	QC_Tools/QCTools_Build/src/IO/ARGUMENTS 
	QC_Tools/QCTools_Build/src/IO/ARGUMENTS/PROPERTIES 
	QC_Tools/QCTools_Build/src/IO/FILE_READERS 
	QC_Tools/QCTools_Build/src/MATRIX 
	QC_Tools/QCTools_Build/src/PARAMETERS 
	QC_Tools/QCTools_Build/src/QC_FUNCTIONS 
	QC_Tools/QCTools_Build/src/STRING_SUPPORT 
	Approx_MVBB/ApproxMVBB/ApproxMVBB_Build/lib
	)

# Generate list of source file
file(
	GLOB_RECURSE
	source_files
	src/*
	)

# Exec declaration
add_executable(
	LOE-CTP-FRAG
	${source_files}
	)

# Add specific build flags for LOE-CTP-FRAG
set_target_properties(
	LOE-CTP-FRAG
	PROPERTIES COMPILE_FLAGS ${BUILD_FLAGS}
	)
	
# Link the libraries
target_link_libraries(
	LOE-CTP-FRAG
	calculateTI 
	IO_SRC 
	ARGUMENTS_SRC 
	PROPERTIES_SRC
	FILE_READERS_SRC
	MATRIX_SRC
 	PARAMETERS_SRC 
	QC_FUNCTIONS_SRC
	STRING_SUPPORT_SRC
 	ApproxMVBB
	)

# Create the script file Load_LOE_CTP_FRAG.sh, useful to set path to use the software
file(
	WRITE 
	scripts/Load_LOE_CTP_FRAG.sh 
"
module load buildenv-intel/2018b-eb
module load vmd/1.9.3-nsc1 
module load Python/3.6.4-nsc2-intel-2018a-eb

export PATH=$PATH:${CMAKE_CURRENT_SOURCE_DIR}/bin:${CMAKE_CURRENT_SOURCE_DIR}/python:${CMAKE_CURRENT_SOURCE_DIR}/scripts
" 
	)
