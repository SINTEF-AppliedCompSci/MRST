# Highlights of MRST 2025a

The new release comes with some improvements in functionality and minor bug fixes. All code from this release can also be cloned or downloaded from the MRST github repository: 

[https://github.com/SINTEF-AppliedCompSci/MRST](https://github.com/SINTEF-AppliedCompSci/MRST)

In particular we would like to highlight the following: 

- New co2-foam modelling for modelling co-injection of CO2 and foam
- New h2store module for the generating PVT data and simulating H2 storage in saline aquifers. 
- New ecpa module implementing the CPA EoS for CCUS fluids to accurately predict phase equilibria of H2O-CO2-H2S-N2-O2-Ar-SO2-CH4-C2H6-C3H8 systems.
- Three phase compositional extension to dual-porosity and dual-porosity-permeability modules.

## New modules

### co2-foam

The co2-foam module provides simulation of CO2 injection into aquifers with mobility control using surfactant. The surfactant can dissolve in both brine and CO2 phases, according to a partitioning coefficient, and adsorbed to reservoir rock. Foam is assumed to be generated wherever the surfactant concentration is large enough. Foam modifies the CO2 mobility according to surfactant concentration, brine saturation and (optionally) the gas velocity.

### h2store

This module provides simulation tools for modeling a hydrogen (H₂) and brine mixture within the Blackoil simulator in MATLAB Reservoir Simulation Toolbox (MRST). It includes tools for implementing the Redlich-Kwong (RK) equation of state (EoS) and generating tabulated PVT data for H₂-brine systems. Additionally, it provides solubility tables derived from ePC-Saft and Henry’s law EoS for precise phase behavior calculations in hydrogen storage simulations. The following features are included: 

- An implementation of the Redlich-Kwong (RK) EoS for H₂-brine mixtures 
- Functionality for precomputing and tabulating PVT Data for Blackoil Simulations 
- Solubility Tables Using ePC-Saft (data) and Henry-Setschnow correlation to get solubility 
- Correlations for water and H2 density calculations 

### ecpa

An new CPA EoS for CCUS fluids to accurately predict phase equilibria of H2O-CO2-H2S-N2-O2-Ar-SO2-CH4-C2H6-C3H8 systems. Contributed by Wei Xiong, Southwest Petroleum University, Chengdu.

A new method is proposed for solving explicit cross-association equations, and the overall CPU time decreases by 70 % for flash calculations and compositional simulations. Extended applicability range of CPA model to CCUS mixtures using vdW mixing rule to replace HV mixing rule. 

References:  

Xiong, Wei, et al. "Phase equilibrium modeling for CCUS fluids using a modified association equation of state." The Journal of Supercritical Fluids (2025): 106543. 

https://doi.org/10.1016/j.supflu.2025.106543 

## Changes to existing modules 

### ad-props
- Added Support for Reading the GSF and WSF Keywords from opm style input.

### dual-porosity
- New compositional functionality added for modeling compositional dual-porosity systems. Contributed by Wei Xiong, Southwest Petroleum University, Chengdu.

### dual-porosity-permeability
- New compositional functionality added for modeling compositional dual-porosity dual-permeability systems. Contributed by Wei Xiong, Southwest Petroleum University, Chengdu.

### test-suite
- Added setup function for the Drogon model.

