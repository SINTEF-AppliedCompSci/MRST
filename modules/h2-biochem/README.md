# MRST Bio-Chemistry Module for Hydrogen Storage Simulation

A comprehensive MATLAB Reservoir Simulation Toolbox (MRST) module for simulating hydrogen storage in depleted reservoirs with bio-chemical reactions and compositional modeling.

## Overview

This module extends MRST's capabilities by integrating a bio-chemistry model with the compositional simulator, specifically designed for hydrogen storage applications. It implements the Soreide-Whitson (SW) equation of state fitted to experimental data and enables simulation of microbial activity affecting hydrogen storage operations.

### Key Features

- **Compositional Modeling**: Full compositional simulation for H₂O-H₂-CO₂-CH₄-N₂ mixtures
- **Bio-Chemistry Integration**: Microbial growth and methanogenesis reactions
- **Equation of State**: Soreide-Whitson EoS fitted to experimental data
- **Bacterial Effects**: Simulation of hydrogen loss due to microbial activity

### Biochemical Reaction Model

The module simulates the methanogenesis reaction:
\[
4\text{H}_2 + \text{CO}_2 \longrightarrow \text{CH}_4 + 2\text{H}_2\text{O} + \text{energy}
\]

For detailed methodology and validation, see our publication:
[**Numerical Modeling of Bio-Reactive Transport During Underground Hydrogen Storage**](https://www.sciencedirect.com/science/article/pii/S0360319925039473)

## Installation

### Prerequisites

- **MATLAB**: Version R2021a or newer
- **MRST**: MATLAB Reservoir Simulation Toolbox (2023b or newer)
- **Required MRST Modules**:
  - `compositional`
  - `ad-blackoil` 
  - `ad-core`
  - `ad-props`
  - `h2store`
  -`biochemistry`
