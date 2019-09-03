# -- Catch tests compiled into target "runtest"


# Copy h5 files to binary dir
FILE(GLOB all_h5_files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.h5)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/${all_h5_files} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})


add_executable(runtest
  bubble.cpp
  runtest.cpp
  function.cpp
  downfold.cpp
  equilibrium.cpp
  integration.cpp
  linalg.cpp
  herm_matrix_member.cpp  
  convolution.cpp
  utilities.cpp
  matsubara.cpp
  dyson.cpp
  test_differentiation.cpp
  )

if (hdf5)
  add_executable(runtest
    fkm_bethe_quench.cpp
  )
endif (hdf5)


if (mpi)
   add_executable(runtest_mpi
	 runtest_mpi.cpp
  	)
endif (mpi)


#include_directories(
#  SYSTEM
#  ../cntr/
#  ../utils/
#  ../hdf5/
#  ../utils/
#  )

# -- Fixme, put this into the Catch framework


target_link_libraries(runtest cntr ${COMMON_LIBS})

if (mpi)
   target_link_libraries(runtest_mpi cntr ${COMMON_LIBS})
endif (mpi)

# -- "test_all" target for running all the test using Catch directly

add_custom_target(test COMMAND runtest -d yes DEPENDS runtest)

if (mpi)
   add_custom_target(test_mpi COMMAND ${MPIEXEC} -np 2 runtest_mpi -d yes DEPENDS runtest_mpi)
endif (mpi)

# ~~ Temp programs

#add_executable(fkm_bethe_quench.ex fkm_bethe_quench.cpp)
#target_link_libraries(fkm_bethe_quench.ex cntr ${COMMON_LIBS})
#install(TARGETS fkm_bethe_quench.ex DESTINATION bin)
