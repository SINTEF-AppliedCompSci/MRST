# computeChargeBalance

Computes the aqueous charge balance of the chemical system.

## SYNOPSIS
~~~~
[state] = chem.computeSurfaceCharges(state)
~~~~

adds the aqueous charge balance  to state as a percentage of total 
concentration of charged species. If 'chargeBalance' is not specified as 
an optional parameter to initState then the value of state.chargeBalance will
likely not be zero.

## OUTPUT
 
### state

state.chargeBalance includes the residual of charge balance as a percentage
of total charge.

## EXAMPLE

~~~~
state = chem.computeSurfaceCharges(state);
chargeRes = chem.getProp(state, 'chargeBalacne');
~~~~

## REQUIRED PARAMETERS
   
### state        

The state structure produced by chem.initstate. Must have the field "components." 