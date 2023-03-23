This module facilitates running ensemble simulations in MRST. 
An ensemble is defined by a 

* baseExample/setup: struct with model, state0, and schedule that defines everything that is common between all ensemble members,
* Stochastic samples: Realizations of uncertain parameters that are unique for each ensemble member
* Quantity of interest (QoI): Description of the simulation output that we aim to estimate the uncertainty of

The module facilitates ensembles that can be used for a wide range of applications,
and has support for (multi-level) Monte Carlo simulations and ensemble-based model calibration.
