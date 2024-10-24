# Simulation tools for H2 storage in aquifers (H2-brine mixture) using Blackoil MRST-Simulator and tabulated PVT data
This module provides simulation tools for modeling a hydrogen (H₂) and brine mixture within the **Blackoil** simulator in **MATLAB Reservoir Simulation Toolbox (MRST)**. 
Our module includes tools for implementing the **Redlich-Kwong** (**RK**) equation of state (EoS) and generating tabulated PVT data for H₂-brine systems. 
Additionally, it provides solubility tables derived from **ePC-Saft** and **Henry’s law** EoS for precise phase behavior calculations in hydrogen 
storage simulations.

## Overview

This module is developed to simulate the phase behavior and thermodynamic properties of hydrogen stored in saline aquifers, with specific consideration of temperature, pressure, and salinity effects on hydrogen solubility and fluid properties.

## Features

- **Implementation of RK Equation of State**: The module includes an implementation of the Redlich-Kwong (RK) EoS for H₂-brine mixtures, allowing for quick thermodynamic calculations that capture real gas behavior and temperature effects.

- ** Scripts for tabulating PVT Data for Blackoil Simulations**: Precomputed PVT tables facilitate efficient H₂-brine mixture simulations in the MRST Blackoil simulator, allowing direct integration with Blackoil simulation workflows.

- **Solubility Tables Using ePC-Saft (data) and Henry-Setschnow correlation**: both can be used to get solubility
- ** Correlations for water and H2 density calculations

## Prerequisites

- **MATLAB R2021a** or newer.
- **MRST** (MATLAB Reservoir Simulation Toolbox) with Blackoil module.

## Installation

1. Clone or download this repository.
2. Add the directory to your MATLAB path:
   ```matlab
   addpath('path/to/H2STORE-module');
## Funding
- This work was funded by Norwegian Research and Innovation Centre for Hydrogen and Ammonia HYDROGENi (Grant No. 333118) and
the Research Council of Norway HYSTORM project (Grant No. 315804)
