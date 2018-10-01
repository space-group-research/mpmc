project(mpmc)
cmake_minimum_required(VERSION 2.8)

option(MPI "Use MPI to parallelize the calculations (requires MPI)" OFF)
option(CUDA "Use CUDA to offload polarization calculations to a GPU (requires CUDA)" OFF)
option(OPENCL "Use OpenCL to offload polarization calculations to a GPU (requires OpenCL)" OFF)
option(QM_ROTATION "Enable Quantum Mechanics Rigid Rotator calculations (requires LAPACK)" OFF)
option(VDW "Enable Coupled-Dipole Van der Waals (requires LAPACK)" OFF)

execute_process(COMMAND bash "-c" "git rev-list HEAD| wc -l |sed 's: ::g'" VERBATIM OUTPUT_VARIABLE REV)
add_definitions(-DVERSION=${REV})

configure_file (
"${PROJECT_SOURCE_DIR}/src/include/cmake_config.h.in"
"${PROJECT_BINARY_DIR}/src/include/cmake_config.h"
)

message("COMPILING FOR ${CMAKE_BUILD_TYPE}")

set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Wall")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Wall")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wall")

#use MATCHES here to use these flags for both regular Clang and AppleClang
if ("${CMAKE_C_COMPILER_ID}" MATCHES "Clang" AND "${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
  #nothing to do here
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU" AND "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  SET( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -static-libstdc++" )
  #elseif ("${CMAKE_COMPILER_ID}" STREQUAL "Intel")
endif()


set(LIB m)

set(INCLUDE src/include src/mersenne ${PROJECT_BINARY_DIR}/src/include)

set(SRC
src/mc/mc_moves.c
src/mc/surface_fit_arbitrary.c
src/mc/qshift.c
src/mc/single_point.c
src/mc/mc.c
src/mc/replay.c
src/mc/pimc.c
src/mc/surface.c
src/mc/surf_fit.c
src/mc/fugacity.cpp
src/mc/cavity.c
src/mc/checkpoint.c
src/histogram/histogram.c
src/energy/lj_buffered_14_7.c
src/energy/bessel.c
src/energy/dreiding.c
src/energy/energy.c
src/energy/polar.c
src/energy/pbc.c
src/energy/disp_expansion.c
src/energy/vdw.c
src/energy/pairs.c
src/energy/bond.c
src/energy/coulombic_gwp.c
src/energy/exp_repulsion.c
src/energy/coulombic.c
src/energy/sg.c
src/energy/lj.c
src/energy/axilrod_teller.cpp
src/mersenne/mt.c
src/mersenne/dSFMT.h
src/mersenne/dSFMT.c
src/mersenne/dSFMT-common.h
src/mersenne/dSFMT-params.h
src/main/quaternion.c
src/main/memnullcheck.c
src/main/main.c
src/main/cleanup.c
src/main/usefulmath.c
src/io/dxwrite.c
src/io/simulation_box.c
src/io/average.c
src/io/output.c
src/io/check_input.c
src/io/input.c
src/io/mpi.c
src/io/read_pqr.c
src/io/setup_ocl.c
src/polarization/thole_field.c
src/polarization/polar_wolf_lookup.c
src/polarization/thole_polarizability.c
src/polarization/thole_matrix.c
src/polarization/polar_ewald.c
src/polarization/thole_iterative.c
)

if(MPI)
	message("-- MPI Enabled")
	find_package(MPI REQUIRED)
	if(NOT MPI_C_FOUND)
		message(FATAL_ERROR "-- MPI not found! Exiting ...")
	endif()
	set(INCLUDE ${INCLUDE} ${MPI_C_INCLUDE_PATH})
	set(LIB ${LIB} ${MPI_C_LIBRARIES})
else()
	message("-- MPI Disabled")
endif()

if(CUDA)
	message("-- CUDA Enabled")
	find_package(CUDA REQUIRED)
	set(SRC ${SRC} src/polarization_gpu/polar_cuda.cu)
else()
	message("-- CUDA Disabled")
endif()

if(OPENCL)
	message("-- OpenCl Enabled")
	set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)
	find_package(OpenCL REQUIRED)
	if(NOT OPENCL_FOUND)
		message(FATAL_ERROR "--OpenCl not found! Exiting ...")
	endif()
	set(SRC ${SRC}
		src/polarization_gpu/polar_ocl.c
		src/io/setup_ocl.c)
	set(INCLUDE ${INCLUDE} ${OpenCL_INCLUDE_DIRS})
	set(LIB ${LIB} ${OpenCL_LIBRARIES})
else()
	message("-- OpenCl Disabled")
endif()

if(QM_ROTATION)
	message("-- QM Rotation Enabled")
	set(SRC ${SRC}
		src/quantum_rotation/rotational_basis.c
		src/quantum_rotation/rotational_eigenspectrum.c
		src/quantum_rotation/rotational_integrate.c
		src/quantum_rotation/rotational_potential.c)
	set(LIB ${LIB} lapack)
else()
	message("-- QM Rotation Disabled")
endif()

if(VDW)
	message("-- CDVDW Enabled")
	if(NOT QM_ROTATION)
		set(LIB ${LIB} lapack)
	endif()
else()
	message("-- CDVDW Disabled")
endif()

include_directories(${INCLUDE})
if(CUDA)
	cuda_add_executable(${PROJECT_NAME} ${SRC})
else()
	add_executable(${PROJECT_NAME} ${SRC})
endif()
target_link_libraries(${PROJECT_NAME} ${LIB})

if(MPI)
	if(MPI_C_COMPILE_FLAGS)
	  set_target_properties(${PROJECT_NAME} PROPERTIES
		COMPILE_FLAGS "${MPI_C_COMPILE_FLAGS}")
	endif()

	if(MPI_C_LINK_FLAGS)
	  set_target_properties(${PROJECT_NAME} PROPERTIES
		LINK_FLAGS "${MPI_C_LINK_FLAGS}")
	endif()
endif()
