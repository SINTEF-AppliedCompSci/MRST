# Simulation Tools for H\(_2\) Storage in Depleted Reservoirs Using Fully Compositional MRST-Simulator and Bio-Chemistry Model

This module integrates a bio-chemistry model with the compositional simulator in the MATLAB Reservoir Simulation Toolbox (MRST). It implements the Soerinde Whitson (SW) equation of state (EoS), which has been fitted to experimental data. The EoS can also generate tabulated PVT data for black-oil simulators.

## Overview

The purpose of this module is to simulate the phase behavior and thermodynamic properties of hydrogen storage in depleted reservoirs. It specifically accounts for the effects of temperature, pressure, and salinity on bacterial activity and hydrogen loss. 

Key features include:
- **SW EoS**: Fitted to the \( \text{H}_2\text{O-H}_2-\text{CO}_2-\text{CH}_4-\text{N}_2 \) mixture for modeling methanogenesis reactions.
- **Biochemical Reaction**: The following reaction is simulated:
  \[
  4\text{H}_2 + \text{CO}_2 \longrightarrow \text{CH}_4 + 2\text{H}_2\text{O} + \text{energy}
  \]

## Prerequisites

To use this module, ensure you have:
- **MATLAB**: Version R2021a or newer.
- **MRST**: MATLAB Reservoir Simulation Toolbox with the compositional module installed.
- **H2STORE**: [Accessible module for hydrogen storage simulation](https://github.com/ElyesAhmed/MRST/tree/hydrogen/modules/H2store).

## Installation

1. Clone or download this repository:
   ```bash
   git clone https://github.com/your-repo-name/bio-chemistry-module.git

# Funding

This work was supported by:

- **HYDROGENi** Norwegian Research and Innovation Centre for Hydrogen and Ammonia (Grant No. 333118).
- **HYSTORM**   Research Council of Norway HYSTORM Project (Grant No. 315804).
