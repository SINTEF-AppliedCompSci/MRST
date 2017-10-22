# computeActivities 

Computes the activities of aqueous species from state using the extended Davies equation. 

## SYNOPSIS
~~~~
[state] = chem.computeActivities(state)
~~~~

adds the field "activities" to state. Activity values can be retrieved using
chem.getProps with 'a' prepended to the species name.

## REQUIRED PARAMETERS
   
### state        

The state structure produced by chem.initstate. Must have the field "components." 

## OUTPUT
 
### state

state.activities includes the activities of aqueous species in units of mol/meter^3.

## EXAMPLE

~~~~
state = chem.computeActivities(state);
aH2O = chem.getProp(state, 'aH2O');
~~~~