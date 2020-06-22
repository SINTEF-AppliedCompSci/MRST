# README #

This README would normally document whatever steps are necessary to get your application up and running.

### What is this repository for? ###

* This is a repo for building and calibrating data driven models for waterflooding optimization.


### How do I get set up? ###

* ADD /ddmodel and /ddmodel/diagnosticsWellPairs in to the PATH in Matlab/MRST

* The examples folder has has Three examples

* These examples follows steps:
    * Settup and simulate a full resolution resorvior model
    * Build a network of flowpath base of flow diagnostics
    * Create the create the ddmodel 
    * Settup parameter and scaling values for the control variables
    * Run initial model
    * Optimize

    
### Overview  of the main functions 
    * createDDmodel is the function that creates the dd model and there are two versions of it:
        createDDmodel_1 creates models with 1 flowpath for interwell conections and createDDmodel can creates moedls qith 2 or more flowpaths. This implementation would be later improved
        
    * Simulate_BFGS is the funciton that evaluates a model, calculates the mistmath between data and prediction, and calculates the gradient of the objective function.
    
    * Parameter in the model has to be scale to the control variables in the optimzer in the interval [0,1]. To help in this scale conversion we used the following functions:
        * problem2control
        * control2problem
        * val2control 
        * control2val 

