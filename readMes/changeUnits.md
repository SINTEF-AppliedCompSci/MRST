# changeUnits

Changes the units of fields in the state variable

## SYNOPSIS

~~~~
[state] = changeUnits(state, fields, units)
~~~~

changes the unit of the the specified field in state. All variables in state are by defualt in SI units. 


## REQUIRED PARAMETERS
   
### state        

The state structure produced by ChemicalModel/initstate. State must be populated by the named variable before it can be retrived.

### field

Name of the field for which the unit change is to occur. Can be a string or cell array of strings.

### unit

Numeric value of unit conversion. If numel(unit) == 1 and field has multple entries then the one value of unit will be applied to all fields.
Otherwise if field is a cell with multple entries then field and unit must have the same size, numel(unit) == numel(field).

## OUTPUT
 
### state

a field with the values of state.field changed to the specified units

## EXAMPLE

change the units of a single field

~~~~
[state] = changeUnits(state, 'activities', mol/litre)
~~~~

change the units of a multiple fields

~~~~
fields = {'acitivities', 'surfaceCharges'}
units = [mol/litre, mili*Coulumb/(nano*meter)^2];

state = chem.getProps(state, fields, units);
~~~~

If state is a cell array of structures the function will loop over each cell. 

The numeric value of the unit conversion has to be defined. The is a bank of unit conversions within the release. 