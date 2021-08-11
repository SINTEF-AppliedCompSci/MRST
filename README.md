# README #

### What is this repository for? ###

* This is a repo for building network models and calibrating them for waterflooding optimization.


### How do I get set up? ###

* ADD /network-model to the PATH in MATLAB/MRST, or just run startup_network_models.m

* The examples folder has has three examples

* These examples follow the same steps:
    * Set up and simulate a full-resolution resorvior model
    * Build a network (i.e., an instance of the Network class)
    * Create the network model (i.e., an instance of the NetworkModel class)
    * Setup parameter and scaling values for the control variables
    * Run initial model
    * Calibrate

    
### Overview  of the main classes 
    * Network: here we create the network and defined the conections between nodes
    * NetworkModel: here we create the MRST model asociated to a given network
