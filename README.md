# README #

### What is this repository for? ###

* This is a repo for building  network models and calibrating it for waterflooding optimization.


### How do I get set up? ###

* ADD /network-model to the PATH in Matlab/MRST, or just run startup_network_models.m

* The examples folder has has Three examples

* These examples follows steps:
    * Settup and simulate a full resolution resorvior model
    * Build a network
    * Create the create the network-model 
    * Settup parameter and scaling values for the control variables
    * Run initial model
    * Calibrate

    
### Overview  of the main classes 
    * Network: here we create the network and defined the conections between nodes
    * NetworkModel: here we create the MRST model asociated to a given network
