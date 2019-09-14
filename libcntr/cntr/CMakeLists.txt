# ~~ The actual target library ~~
add_library(cntr SHARED
    fourier.cpp
    linalg_eigen.cpp
    integration.cpp
    integration_extern_templates.cpp
    cntr_convolution_extern_templates.cpp
    cntr_pseudo_convolution_extern_templates.cpp
    #cntr_cyclic_timestep_extern_templates.cpp
    cntr_differentiation_extern_templates.cpp
    cntr_dyson_extern_templates.cpp
    cntr_dyson_omp_extern_templates.cpp
    cntr_equilibrium_extern_templates.cpp
    #cntr_full_matrix_extern_templates.cpp
    cntr_function_extern_templates.cpp
    cntr_herm_matrix_extern_templates.cpp
    cntr_herm_matrix_timestep_extern_templates.cpp
    cntr_herm_matrix_timestep_view_extern_templates.cpp
    cntr_herm_pseudo_extern_templates.cpp
    cntr_equilibrium_extern_templates.cpp
    cntr_utilities_extern_templates.cpp
    cntr_vie2_extern_templates.cpp
    cntr_distributed_array_extern_templates.cpp
    cntr_distributed_timestep_array_extern_templates.cpp
    cntr_getset_extern_templates.cpp
    )
set_property(TARGET cntr APPEND_STRING PROPERTY COMPILE_FLAGS " -w")

if (hdf5)
    add_subdirectory(hdf5)
    target_link_libraries(cntr cntr_hdf5)
endif (hdf5)

# Add all headers as public-headers to `libcntr`. That way they will
# automatically be installed alongside `libcntr`.
file(GLOB cntr_HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
set_target_properties(cntr PROPERTIES PUBLIC_HEADER "${cntr_HEADERS}")

# ~~ Installation targets ~~
install(TARGETS cntr
        LIBRARY DESTINATION lib
        PUBLIC_HEADER DESTINATION include/cntr)
if (hdf5)
    install(DIRECTORY hdf5 DESTINATION include/cntr
            FILES_MATCHING PATTERN "*.hpp")
endif (hdf5)
install(DIRECTORY cntr/utils DESTINATION include/cntr
        FILES_MATCHING PATTERN "*.hpp")
install(FILES utils/read_inputfile.hpp DESTINATION include/cntr/utils)	


