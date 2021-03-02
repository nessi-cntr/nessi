set( regular_EXECUTABLES
     test_hodlr.x
)

set( EXE_test_hodlr.x_SOURCES test_hodlr.cpp )

# ~~ Paths and Subdirs ~~
include_directories(${CMAKE_INCLUDE_PATH})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# ~~ Libraries ~~
link_directories(${CMAKE_LIBRARY_PATH})

foreach( exe ${regular_EXECUTABLES} )
    # declare the target...
    add_executable( ${exe} ${EXE_${exe}_SOURCES} )
    target_link_libraries(${exe} cntr ${COMMON_LIBS})
endforeach()

if(hdf5)
    foreach( exe ${withhdf5_EXECUTABLES} )
        # declare the target...
        add_executable( ${exe} ${EXE_${exe}_SOURCES} )
        target_link_libraries(${exe} cntr cntr_hdf5 ${COMMON_LIBS})
    endforeach()
endif(hdf5)

if(mpi AND hdf5) 
    foreach( exe ${mpi_EXECUTABLES} )
        # declare the target...
        add_executable( ${exe} ${EXE_${exe}_SOURCES} )
        target_link_libraries(${exe} cntr cntr_hdf5 ${COMMON_LIBS})
    endforeach()
endif(mpi AND hdf5)