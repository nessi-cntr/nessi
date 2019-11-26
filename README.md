[![Build Status](https://travis-ci.org/nessi-cntr/nessi.svg?branch=master)](https://travis-ci.org/nessi-cntr/nessi/)
# Welcome to the Non-Equilibrium Systems Simulation (NESSi) package!

## What is NESSi?

**NESSi** is an open-source software package for the manipulation of nonequilibrium Green’s functions defined on the Kadanoff-Baym contour. The Green's function method in its time-dependent formulation is a versatile framework for the solution of interacting many-body problems out of equilibrium. **NESSi** provides classes representing the various types of Green’s functions, implements the basic operations on these functions and allows to solve the corresponding equations of motion. The library is aimed at the study of transient dynamics from an initial equilibrium state, induced by time-dependent model parameters.

**Overview:**
* NESSi provides tools for constructing Feynman diagram and solving equations of motion for non-equilibrium Green's functions on the Kadanoff-Baym contour
* NESSi is based on high-order quadrature rules: for *N* time slices, the error scales like *O(N<sup>-p</sup>)* with *p* up to 7.
* Efficient distributed-memory parallelization over reciprocal space allows large-scale calculations on extended systems.

## Content

The NESSi program package contains two major directories:

Directory | Content
------------ | -------------
[libcntr](libcntr/) | The **libcntr** library.
[examples](examples/) | Example programs based on **libcntr**.

More information on the directory structure and the source files
can be found in the respective directories.

## Documentation and tutorials

The full documentation, installation instructions, examples and tutorials is available [here](http://www.nessi.tuxfamily.org).

## How to cite

Please cite the following paper whenever you use parts of **NESSi**: 

Michael Schüler, Denis Golež, Yuta Murakami, Nikolaj Bittner, Andreas Herrmann, Hugo U. R. Strand, Philipp Werner, Martin Eckstein, [arXiv:1911.01211 [cs]](http://arxiv.org/abs/1911.01211)

## License

This project is licensed under the MPL 2.0 License - see the [LICENSE.md](LICENSE.md) file for details
