# A list of executables to build

set( regular_EXECUTABLES
     hubbard_chain_2b.x
     hubbard_chain_gw.x
     hubbard_chain_tpp.x
     test_equilibrium.x
     test_nonequilibrium.x
     integration.x
     test_jason.x
)

set( withhdf5_EXECUTABLES
     Holstein_bethe_Nambu_Migdal.x
     Holstein_bethe_Nambu_uMig.x
     Holstein_bethe_Migdal.x
     Holstein_bethe_uMig.x
     Holstein_impurity_singlebath_Migdal.x
     Holstein_impurity_singlebath_uMig.x
     trunc_bethe_start.x
     trunc_bethe.x
)

set (mpi_EXECUTABLES
    gw.x
)

# Add NESS-related executable if ness = ON; all require also hdf5=ON
set( ness_EXECUTABLES
     ness2_siam.x
     ness2_bethe.x
     ness2_bethe_prop.x
     ness2_bethe_trunc.x     
)

set( EXE_hubbard_chain_2b.x_SOURCES hubbard_chain_selfen_impl.cpp hubbard_chain_2b.cpp )
set( EXE_hubbard_chain_gw.x_SOURCES hubbard_chain_selfen_impl.cpp hubbard_chain_gw.cpp )
set( EXE_hubbard_chain_tpp.x_SOURCES hubbard_chain_selfen_impl.cpp hubbard_chain_tpp.cpp )


set( EXE_test_equilibrium.x_SOURCES test_equilibrium.cpp )
set( EXE_test_nonequilibrium.x_SOURCES test_nonequilibrium.cpp )
set( EXE_test_jason.x_SOURCES test_jason.cpp )
set( EXE_gw.x_SOURCES gw.cpp gw_latt_impl.cpp gw_kpoints_impl.cpp  gw_selfene_impl.cpp )
set( EXE_integration.x_SOURCES integration.cpp )
set( EXE_Holstein_bethe_Nambu_Migdal.x_SOURCES Holstein_impurity_impl.cpp Holstein_utils_impl.cpp Holstein_bethe_Nambu_Migdal.cpp )
set( EXE_Holstein_bethe_Nambu_uMig.x_SOURCES Holstein_impurity_impl.cpp Holstein_utils_impl.cpp Holstein_bethe_Nambu_uMig.cpp )

set( EXE_Holstein_bethe_Migdal.x_SOURCES Holstein_impurity_impl.cpp Holstein_utils_impl.cpp Holstein_bethe_Migdal.cpp )
set( EXE_Holstein_bethe_uMig.x_SOURCES Holstein_impurity_impl.cpp Holstein_utils_impl.cpp Holstein_bethe_uMig.cpp )
set( EXE_Holstein_impurity_singlebath_Migdal.x_SOURCES Holstein_impurity_impl.cpp Holstein_utils_impl.cpp Holstein_impurity_singlebath_Migdal.cpp )
set( EXE_Holstein_impurity_singlebath_uMig.x_SOURCES Holstein_impurity_impl.cpp Holstein_utils_impl.cpp Holstein_impurity_singlebath_uMig.cpp )


# demo programs for ness2 and trunc:
set( EXE_ness2_bethe_prop.x_SOURCES ness2_bethe_prop.cpp )
set( EXE_ness2_bethe_trunc.x_SOURCES ness2_bethe_trunc.cpp )
set( EXE_ness2_siam.x_SOURCES ness2_siam.cpp )
set( EXE_ness2_bethe.x_SOURCES ness2_bethe.cpp )
set( EXE_trunc_bethe.x_SOURCES trunc_bethe.cpp )
set( EXE_trunc_bethe_start.x_SOURCES trunc_bethe_start.cpp )


# ~~ Paths and Subdirs ~~
include_directories(${CMAKE_INCLUDE_PATH})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# ~~ Libraries ~~
link_directories(${CMAKE_LIBRARY_PATH})
find_library(CNTR_LIB cntr PATHS ${CMAKE_LIBRARY_PATH}
NO_DEFAULT_PATH)
message(STATUS "cntr library found: ${CNTR_LIB}")


# ~~ NESSi 2.0: added libomp in this part (macOS) ~~
foreach( exe ${regular_EXECUTABLES} )
    # declare the target...
    add_executable( ${exe} ${EXE_${exe}_SOURCES} )
    if (omp)
        target_link_libraries(${exe} cntr ${COMMON_LIBS} OpenMP::OpenMP_CXX)
    else(omp)
        target_link_libraries(${exe} cntr ${COMMON_LIBS})
    endif(omp)
endforeach()

if(hdf5)
    foreach( exe ${withhdf5_EXECUTABLES} )
        # declare the target...
        add_executable( ${exe} ${EXE_${exe}_SOURCES} )
        if (omp)
            target_link_libraries(${exe} cntr cntr_hdf5 ${COMMON_LIBS} OpenMP::OpenMP_CXX)
        else(omp)
            target_link_libraries(${exe} cntr cntr_hdf5 ${COMMON_LIBS})
        endif(omp)
    endforeach()
endif(hdf5)

if(mpi AND hdf5) 
    foreach( exe ${mpi_EXECUTABLES} )
        # declare the target...
        add_executable( ${exe} ${EXE_${exe}_SOURCES} )
        if (omp)
            target_link_libraries(${exe} cntr cntr_hdf5 ${COMMON_LIBS} OpenMP::OpenMP_CXX)
        else(omp)
            target_link_libraries(${exe} cntr cntr_hdf5 ${COMMON_LIBS})
        endif(omp)
    endforeach()
endif(mpi AND hdf5)

if (ness AND hdf5)
    foreach( exe ${ness_EXECUTABLES} )
        add_executable( ${exe} ${EXE_${exe}_SOURCES} )
        if (omp)
            target_link_libraries(${exe} cntr cntr_hdf5 ${COMMON_LIBS} OpenMP::OpenMP_CXX)
        else(omp)
            target_link_libraries(${exe} cntr cntr_hdf5 ${COMMON_LIBS})
        endif(omp)
endforeach()
endif()

