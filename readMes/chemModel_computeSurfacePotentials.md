# computeSurfacePotentials

Computes the surface potential of each layer of each surface in the chemical system.

## SYNOPSIS
~~~~
[state] = chem.computeSurfacePotentials(state)
~~~~

adds the potential of each layer of each surface to the field "surfacePotentials" in state. 
Surface potential values can be retrieved using the getProp command with the
surface functional group, followed by '_Psi_' followed by the layer number.

## OUTPUT
 
### state

state.surfaceCharges includes the charge of each surface and layer in Volts.

## EXAMPLE

~~~~
state = chem.computeSurfacePotentials(state);
potential0 = chem.getProp(state, '>SiO_Psi_0');
potential1 = chem.getProp(state, '>SiO_Psi_1');
~~~~

## REQUIRED PARAMETERS
   
### state        

The state structure produced by chem.initstate. Must have the field "components." 