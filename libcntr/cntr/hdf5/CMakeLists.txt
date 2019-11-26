
# ~~ cntr_hdf5 shared library 

set(cntr_hdf5_src hdf5_interface.cpp)
add_library(cntr_hdf5 SHARED ${cntr_hdf5_src})
#target_link_libraries(cntr_hdf5 cntr ${COMMON_LIBS})
target_link_libraries(cntr_hdf5 ${COMMON_LIBS} ${HDF5_LIB})
install(TARGETS cntr_hdf5 DESTINATION lib)

