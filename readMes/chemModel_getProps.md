# getProp/s 

Allows the retreval of values from the state variables created by the solution of the chemical system.

## SYNOPSIS
~~~~
[value] = chem.getProp(state, name)
~~~~

retrieves the value of the named property from the state variable

~~~~
[values1, value2] = chem.getProps(state, name1, name2)
~~~~

retrieves values of multiple properties from state

## OUTPUT
 
### value

vector value of named variable within state

## EXAMPLE

Grab the concentration of H+ in the system

~~~~
H = chem.getProp(state, 'H+');
~~~~

Grab the concentration and acitivity of H+

~~~~
[H, aH] = chem.getProps(state, 'H+','aH+');
~~~~



## REQUIRED PARAMETERS
   
### state        

The state structure produced by ChemicalModel/initstate. State must be populated by the named variable before it can be retrived.

### name

Name of the vairble to retrieve, string. 
