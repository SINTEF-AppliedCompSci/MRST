The core module offers data structures and routines for creating and manipulating:
  - geological description: structured and unstructured grids

  - petrophysical properties (porosity, permeability, net to gross, etc)

  - drive mechanisms (wells, boundary conditions)

  - reservoir state (pressures, saturations, fluxes, etc.)

This includes, in particular, a large number of grid-factory routines and a routine for computing transmissibilities.


In addition, the core module provides:

  - a library for automatic differentiation geared towards sparse matrices and coupled PDEs

  - plotting of cell and face data defined over MRST grids

  - physical units (length, time, mass, pressures, etc) and conversion routines

  - basic reading, parsing, and writing ECLIPSE input data

  - various utility routines and functionality


To get started we suggest you look at some examples. To list all the examples in mrst core at the matlab command prompt type:

>> mrstExamples

Each module also has its own examples which can also be listed, e.g. to list examples for the ad-core module you could type:

>> mrstExamples ad-core




