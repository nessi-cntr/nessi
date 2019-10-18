# Welcome to the Non-Equilibrium Systems Simulation (NESSi) package!

## What is NESSi?

*NESSi* is an open-source software package for the manipulation of nonequilibrium Green’s functions defined on the Kadanoff-Baym contour. The Green's function method in its time-dependent formulation is a versatile framework for the solution of interacting many-body problems out of equilibrium. *NESSi* provides classes representing the various types of Green’s functions, implements the basic operations on these functions and allows to solve the corresponding equations of motion. The library is aimed at the study of transient dynamics from an initial equilibrium state, induced by time-dependent model parameters.

**Overview:**
* NESSi provides tools for constructing Feynman diagram and solving equations of motion for non-equilibrium Green's functions on the Kadanoff-Baym contour
* NESSi is based on high-order quadrature rules: for *N* time slices, the error scales like *O(N<sup>-p</sup>)* with *p* up to 7.
* Efficient distributed-memory parallelization over reciprocal space allows large-scale calculations on extended systems.

## Documentation and tutorials

The full documentation, installation instructions, examples and tutorials is available [here](http://www.nessi.tuxfamily.org).

## How to cite

Please cite the following paper whenever you use parts of *NESSi*: 

Michael Schüler, Denis Golež, Yuta Murakami, Nikolaj Bittner, Andreas Herrmann, Hugo U. R. Strand, Philipp Werner, Martin Eckstein, CPC XX, XX (2019)

This project is licensed under the MPL License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to anyone whose code was used
* Inspiration
* etc

