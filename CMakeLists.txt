project(mpmc)
cmake_minimum_required(VERSION 2.8)

#if (CMAKE_VERSION VERSION_LESS "3.1")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -D_GLIBCXX_USE_CXX11_ABI=0")
#else ()
  #set (CMAKE_CXX_STANDARD 11)
#endif()


option(MPI "Use MPI to parallelize the calculations (requires MPI)" OFF)
option(CUDA "Use CUDA to offload polarization calculations to a GPU (requires CUDA)" OFF)
option(QM_ROTATION "Enable Quantum Mechanics Rigid Rotator calculations (requires LAPACK)" OFF)
option(VDW "Enable Coupled-Dipole Van der Waals (requires LAPACK)" OFF)

execute_process(COMMAND bash "-c" "git rev-list HEAD| wc -l |sed 's: ::g'" VERBATIM OUTPUT_VARIABLE REV)
add_definitions(-DVERSION=${REV})

configure_file (
"${PROJECT_SOURCE_DIR}/src/include/cmake_config.h.in"
"${PROJECT_BINARY_DIR}/src/include/cmake_config.h"
)

message("Compiling For ${CMAKE_BUILD_TYPE}")

set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -g -Wall")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -g -Wall")

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -g -Wall")



#use MATCHES here to use these flags for both regular Clang and AppleClang
if ("${CMAKE_C_COMPILER_ID}" MATCHES "Clang" AND "${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
  #nothing to do here
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU" AND "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  SET( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -static-libstdc++" )
elseif ("${CMAKE_COMPILER_ID}" STREQUAL "Intel") #just in case
endif()


set(LIB m)

set(INCLUDE src/include ${PROJECT_BINARY_DIR}/src/include)

set(SRC
src/mc/mc_moves.c
src/mc/surface_fit_arbitrary.c
src/mc/surface_multi_fit.c
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
src/io/histogram.c
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
src/main/quaternion.c
src/main/memnullcheck.c
src/main/main.c
src/main/cleanup.c
src/main/usefulmath.c
src/main/rand.c
src/io/dxwrite.c
src/io/simulation_box.c
src/io/average.c
src/io/output.c
src/io/check_input.c
src/io/input.c
src/io/mpi.c
src/io/read_pqr.c
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
    set(SRC ${SRC} src/polarization/polar_cuda_pcg.cu)
    set(CUDA_PROPAGATE_HOST_FLAGS OFF)
    find_package(Threads REQUIRED)
else()
    message("-- CUDA Disabled")
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
    cuda_add_cublas_to_target(${PROJECT_NAME})
    target_link_libraries(${PROJECT_NAME} ${CMAKE_THREAD_LIBS_INIT})
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
