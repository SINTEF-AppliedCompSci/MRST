# computeActivities 

Computes the acitivities of aqueous species from state using the extended Davies equaiton. 

## SYNOPSIS
~~~~
[state] = chem.computeActivities(state)
~~~~

adds the field "activities" to state. Activity values can be retrieved using
chem.getProps with 'a' prepended to the speices name.

## OUTPUT
 
### state

state.activities includes the activites of aqueous species in units of mol/meter^3.

## EXAMPLE

~~~~
state = chem.computeActivities(state);
aH2O = chem.getProp(state, 'aH2O');
~~~~

## REQUIRED PARAMETERS
   
### state        

The state structure produced by chem.initstate. Must have the field "components." 