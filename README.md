# README #

The matlab-geoChemistry repository contains tools for the solution of equilibrium geochemical systems including aqueous and surface chemistry for use in batch and transport settings. 

### Summary ###

The main function of this repository, ChemicalModel.m, allows the creation and solution of arbitrarily complex aqueous chemistry systems including a number of surface chemistry models assuming local chemical equilibrium. The function leverages the tools developed by the SINTEF [MRST team](http://www.sintef.no/projectweb/mrst/) including [automatic differentiation](https://en.wikipedia.org/wiki/Automatic_differentiation). The chemical model created can be used to calculate batch reaction or can be coupled to flow within MRST.

Models supported:

* aqueous speciation
* aqueous activity
* surface chemistry
    * triple layer model
    * diffuse layer model
    * basic stern model
    * constant capacitance model
    * ion exchange

Supported soon:

* dissolution/precipitation

### Installation ###

1. Install [MRST](http://www.sintef.no/projectweb/mrst/downloadable-resources/). 
2. Add the matlab-geoChemistry folder to the module folder of MRST.
3. Create a file named startup_user.m within the MRST folder, at the same level as startup.m.
4. In startup_user.m add the line
~~~~
mrstPath('register', 'geochemistry', 'path/to/repo/matlab-geochemistry')
~~~~

### Use ###

Once MRST is installed and made aware of the location of matlab-geoChemistry the module can be used like any other MRST module. 

Before any script that relies on the repository is run, MRST must be started. This is done by running the file startup.m which is loacted inside of your MRST directory.

To use the geochemistry module in a Matlab script include the command

~~~~~
mrstModule add geochemistry
~~~~~

this will make the contents of the geochemistry directory available in the workspace.

## Functionality ##

[ChemModel.m](readMes/ChemModel.md)
Used for the instatiation of the chemical object, making all embedded functions available. 

* [ChemModel/initState](readMes/chemModel_initState.md) Used to solve the chemical system based on specified user inputs

* [ChemModel/computeActivities](readMes/chemModel_computeActivities.md) adds the acitivities of aqueous species to state

* [ChemModel/computeChargeBalance](readMes/chemModel_computeChargeBalance.md) adds the residual of the charge balance to state

* [ChemModel/computeSurfacePotentials](readMes/chemModel_computeSurfacePotentials.md) adds the potential of each surface and layer to state

* [ChemModel/computeSurfaceCharges](readMes/chemModel_computeSurfaceCharges.md) adds the charge of each surface and layer to state

[mergeChemicalModels.m](readMes/mergeChemicalModels.md) combine two or more chemical models

[changeUnits.m](readMes/changeUnits.md) changes the units of the specified quantity whatever you want