Build instructions
==================

To build the tests eigen must be installed in the include path (under ./eigen3/). The include and library paths can be specified when running cmake as

```
    # From the root directory .../libcntr/
    # Create build directory
    mkdir cbuild
    cd cbuild
```

    # Cmake configuration step (replace paths according to your system)

```
    cmake \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_INCLUDE_PATH=/opt/local/include \
        -DCMAKE_LIBRARY_PATH=/opt/local/lib \
        ..
```

    # Build the tests
```
    make -j
```

    # Executables now reside in ./cbuild/cntr/
    cd ./cntr
    ./test_vie2 1.0 1.0 10 100 5

HDF5
====

To enable hdf5 support use the `hdf5` option in cmake

    cmake -Dhdf5=ON

MPI and OpenMP
==============

To turn on OpenMP and/or MPI parallelization define the options `omp` and/or `mpi` in the cmake step, respectively, and specify your c and c++ mpi compilers according to

    CC=mpicc CXX=mpic++ \
    cmake \
        -DCMAKE_BUILD_TYPE=Release \
        -Domp=ON \
        -Dmpi=ON \
        -DCMAKE_INCLUDE_PATH=/opt/local/include \
        -DCMAKE_LIBRARY_PATH=/opt/local/lib \
        ..

Note, that CMake should detect the correct options in most cases. Hence, it
should not be necessary to explicitely define `CC`, and `CXX`.

doxygen
==============

To generate documentation run

    doxygen Doxyfile.txt

it creates a html directory inside an output doc directory. In html directory there is an index.html, which contains the generated HTML documentation.

To make a description in the code use the following comment block starting with two *'s:

    /** \brief Here is a short documentation.
    * Detailed documentation starts here.
    */

Class documentation starts with \class. The first sentence after \class is a short description of the class.

    /** \class here is a short class documentation.
    * Detailed documentation starts here.
    */


To describe parameters of a function use:

    @param  some_parameter  Description of some parameter

To include formulas in the description use:

    \f{equation*}{
        some_formula_with_latex_commands
    }

Example of documentation:

    /** \brief <b> Left-Mixing convolution at given time-step. </b>
      *
      * <!-- ====== DOCUMENTATION ====== -->
      *
      *   \par Purpose
      * <!-- ========= -->
      *
      * > Here we calculate
      * > \f{eqnarray*}{
      *     C &=& A*B\\
      *     C^{tv}(t,\tau') &=& \int_0^\beta ds  A^{tv}(t,s) B^M(s - \tau')
      *                       + \int_0^t ds A^R(t,s) B^{tv}(s,\tau') \f}
      * >
      * > At time-step `t = n h` for all \f$\tau'\f$.
      * > The result is written into `ctv`.
      * >
      * > `ctv` is a pointer to the timestep of elements of type GG,
      * > so it must have size `(C.ntau_ + 1) * C.element_size_`
      *
      * <!-- ARGUMENTS
      *      ========= -->
      *
      * @param beta
      * > inverse temperature
      *
      * @param h
      * > time step
      */


Tests
=====

To build the tests activate the option `test`. The executables will be placed in `./cbuild/cntr/`. Those tests that depend on OpenMP will only be built if the `omp` option is activated.

    cmake \
        -Dtest=ON
        -Domp=ON \
        -Dmpi=ON \
        -DCMAKE_INCLUDE_PATH=/opt/local/include \
        -DCMAKE_LIBRARY_PATH=/opt/local/lib \
        ..

Wiki
=====

The wiki is a git repository. To clone it, use: 

``` 
git clone https://bitbucket.org/computationalphysicsunifr/libcntr.git/wiki
```
