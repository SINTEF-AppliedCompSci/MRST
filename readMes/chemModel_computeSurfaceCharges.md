# computeSurfaceCharges 

Computes the surface charge of each layer of each surface in the chemical system.

## SYNOPSIS
~~~~
[state] = chem.computeSurfaceCharges(state)
~~~~

adds the charge of each layer of each surface to the field "surfaceCharges" in state. 
Surface charge values can be retrieved using the getProp command with the
surface functional group, followed by '_sig_' followed by the layer number.

## OUTPUT
 
### state

state.surfaceCharges includes the charge of each surface and layer in C/meter^2.

## EXAMPLE

~~~~
state = chem.computeSurfaceCharges(state);
charge0 = chem.getProp(state, '>SiO_sig_0');
charge1 = chem.getProp(state, '>SiO_sig_1');
~~~~

## REQUIRED PARAMETERS
   
### state        

The state structure produced by chem.initstate. Must have the field "components." 