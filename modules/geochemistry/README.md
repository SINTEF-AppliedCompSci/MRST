# README #

The matlab-geoChemistry repository contains tools for the solution of equilibrium geochemical systems including aqueous and surface chemistry for use in batch and transport settings. 

### Summary ###

The main function of this repository, ChemicalModel, allows the creation and solution of arbitrarily complex chemistry systems including equilibrium aqueous speciation, surface complexation, ion exchange, redox chemistry, dissolution/precipitation and equilbrium with gas phases. The function leverages the tools developed by the SINTEF [MRST team](http://www.sintef.no/projectweb/mrst/) including [automatic differentiation](https://en.wikipedia.org/wiki/Automatic_differentiation). The chemical model created can be used to calculate batch reaction or can be coupled to flow within MRST. See the user guide pdf for details on use, equations, examples, and benchmarks. 

Models supported:

* aqueous speciation
* aqueous activity
* surface chemistry
    * triple layer model
    * diffuse layer model
    * basic stern model
    * constant capacitance model
    * ion exchange
* dissolution/precipitation
* gas phase equilibrium

### Installation ###

1. Install [MRST](http://www.mrst.no). 
2. Add the matlab-geoChemistry folder to the "modules" folder of MRST.
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

