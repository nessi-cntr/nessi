# A list of executables to build
set( project_EXECUTABLES
     hubbard_chain_2b.x
     hubbard_chain_gw.x
     hubbard_chain_tpp.x
     test_equilibrium.x
     test_nonequilibrium.x
     integration.x
     Holstein_bethe_Nambu_Migdal.x
     Holstein_bethe_Nambu_uMig.x
     Holstein_bethe_Migdal.x
     Holstein_bethe_uMig.x
     Holstein_impurity_singlebath_Migdal.x
     Holstein_impurity_singlebath_uMig.x
     gw.x
)

set( EXE_hubbard_chain_2b.x_SOURCES hubbard_chain_selfen_impl.cpp hubbard_chain_2b.cpp )
set( EXE_hubbard_chain_gw.x_SOURCES hubbard_chain_selfen_impl.cpp hubbard_chain_gw.cpp )
set( EXE_hubbard_chain_tpp.x_SOURCES hubbard_chain_selfen_impl.cpp hubbard_chain_tpp.cpp )

set( EXE_test_equilibrium.x_SOURCES test_equilibrium.cpp )
set( EXE_test_nonequilibrium.x_SOURCES test_nonequilibrium.cpp )
set( EXE_gw.x_SOURCES gw.cpp gw_latt_impl.cpp gw_kpoints_impl.cpp  gw_selfene_impl.cpp )
set( EXE_integration.x_SOURCES integration.cpp )
set( EXE_Holstein_bethe_Nambu_Migdal.x_SOURCES Holstein_impurity_impl.cpp Holstein_utils_impl.cpp Holstein_bethe_Nambu_Migdal.cpp )
set( EXE_Holstein_bethe_Nambu_uMig.x_SOURCES Holstein_impurity_impl.cpp Holstein_utils_impl.cpp Holstein_bethe_Nambu_uMig.cpp )

set( EXE_Holstein_bethe_Migdal.x_SOURCES Holstein_impurity_impl.cpp Holstein_utils_impl.cpp Holstein_bethe_Migdal.cpp )
set( EXE_Holstein_bethe_uMig.x_SOURCES Holstein_impurity_impl.cpp Holstein_utils_impl.cpp Holstein_bethe_uMig.cpp )
set( EXE_Holstein_impurity_singlebath_Migdal.x_SOURCES Holstein_impurity_impl.cpp Holstein_utils_impl.cpp Holstein_impurity_singlebath_Migdal.cpp )
set( EXE_Holstein_impurity_singlebath_uMig.x_SOURCES Holstein_impurity_impl.cpp Holstein_utils_impl.cpp Holstein_impurity_singlebath_uMig.cpp )


# ~~ Paths and Subdirs ~~
include_directories(${CMAKE_INCLUDE_PATH})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# ~~ Libraries ~~
link_directories(${CMAKE_LIBRARY_PATH})

foreach( exe ${project_EXECUTABLES} )
    # declare the target...
    add_executable( ${exe} ${EXE_${exe}_SOURCES} )
    target_link_libraries(${exe} cntr cntr_hdf5 ${COMMON_LIBS})
endforeach()

