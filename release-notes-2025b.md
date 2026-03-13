# Highlights of MRST 2025b

The new release comes with some improvements in functionality and minor bug fixes. All code from this release can also be cloned or downloaded from the MRST github repository: 

[https://github.com/SINTEF-AppliedCompSci/MRST](https://github.com/SINTEF-AppliedCompSci/MRST)

In particular we would like to highlight the following: 

- A new module for modeling microbial biogeochemical reactions in Underground Hydrogen Storage (UHS).
- Søreide–Whitson Equation of State for H₂ phase behaviour in the compositional module.
- Added functionality for using the hysteresis from co2lab-mit with the TwoPhaseGasWaterModel.
- Improved handling of very low mass flow rates in the WellboreModel.
- Functionality for using the Embedded Discrete Fracture Model (EDFM) with Corner-Point grids.


## New modules

### h2-biochem

A new module for modeling microbial biogeochemical reactions in Underground Hydrogen Storage (UHS). The module provides a fully implicit, fully coupled framework that extends MRST's compositional simulator to account for hydrogen-consuming microbes (e.g., methanogenic archaea).

Key features include:
- (1) Microbial growth/decay via a double Monod formulation 
- (2) Søreide–Whitson Equation of State for accurate H₂ phase behavior 
- (3) Modeling of bio-clogging and salinity effects on H2 solubility.

Contributed by Elyes Ahmed (Sintef Digital), Stéphanie Delage Santacreu (Université de Pau et des Pays de l’Adour, Pau France).
Reference: Elyes Ahmed, Brahim Amaziane, Salaheddine Chabab, Stéphanie Delage Santacreu, Guillaume Galliéro, Olav Møyner, Xavier Raynaud, Modeling and simulation of coupled biochemical and two-phase compositional flow in underground hydrogen storage, International Journal of Hydrogen Energy, Volume 168, 2025, 150947, ISSN 0360-3199, https://doi.org/10.1016/j.ijhydene.2025.150947.

## Changes to existing modules 

### co2lab
- Improved and more precise inventory computation in 2D and 3D (postprocessStates, masssTrappingDistribution).  In particular, should handle the capillary fringe case in a more consistent way.
- Improved convergence when using the hysteresis in co2lab-mit, and added functionality for using the hysteresis from co2lab-mit with the TwoPhaseGasWaterModel in co2lab
- More robust handling of the capillary fringe model
- General bugfixes

### compositional
- Added Søreide–Whitson Equation of State for modelling H₂ phase behaviour

### geothermal
This PR introduces minor updates to the geothermal module to improve handling of very low mass flow rates in the WellboreModel and adds support for CompositeModel instances in TestCase.

- Enhanced WellboreModel flux equations to handle very low friction loss scenarios by adding mass flux velocity term
- Extended TestCase class to properly support CompositeModel instances for plotting and driving force detection
- Improved wellbore trajectory processing and geothermal equilibrium initialization

### hfm
- Functionality for using the Embedded Discrete Fracture Model (EDFM) with Corner-Point grids. Contributed by Mohamad Karimi, Petroleum Reservoir Engineering, Iran University of Science and Technology