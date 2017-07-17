function [ states ] = changeUnits( states, field, unit_conv )
%changeUnits converts the specified field in states to the specified units. 
% 
% SYNOPSIS:
% 
%         [state] = changeUnits(state, field, unit)
% 
% OUTPUT:
%  
%	state           -a field with the values of state.field changed to the specified units
% 
% REQUIRED PARAMETERS:
%
%	state           -The state structure produced by ChemicalModel/initstate. 
%                       State must be populated by the named variable before
%                       it can be retrived.
%
%	field           -Name of the field for which the unit change is to occur. 
%                       Can be a string or cell array of strings.
%
%	unit            -Numeric value of unit conversion. If numel(unit) == 1 
%                       and field has multple entries then the one value of 
%                       unit will be applied to all fields. Otherwise if field
%                       is a cell with multple entries then field and unit must 
%                       have the same size, numel(unit) == numel(field).unit_conv = unit_conv.^-1;
%
% EXAMPLE:
% 
%   change the units of a single field
% 
%         [state] = changeUnits(state, 'activities', mol/litre)
% 
%   change the units of a multiple fields
% 
%         fields = {'acitivities', 'surfaceCharges'}
%         units = [mol/litre, mili*Coulumb/(nano*meter)^2];
% 
%         state = chem.getProps(state, fields, units);
% 
%   If state is a cell array of structures the function will loop over each cell. 
% 
%   The numeric value of the unit conversion must be defined. The is a bank 
%   of unit conversions within the release. 
% 

assert(isnumeric(unit_conv), 'The variable "unit_conv" must be numeric.');
unit_conv = unit_conv.^-1;

if numel(states) == 1
    
    if ~iscell(field)
        states.(field) = states.(field)*unit_conv;
    else
        assert(numel(field) == numel(unit_conv) | numel(unit_conv) == 1, 'The number of elements in unit_conv must be 1 or equal to the number of cells in field.');
        if numel(unit_conv) == numel(field)
            for i = 1 : numel(field)
                states.(field{i}) = states.(field{i})*unit_conv(i);
            end
        else
            for i = 1 : numel(field)
                states.(field{i}) = states.(field{i})*unit_conv;
            end
        end
    end
        
else

    for i = 1 : numel(states)
    
        if ~iscell(field)
            states{i}.(field) = states{i}.(field)*unit_conv;
        else
            assert(numel(field) == numel(unit_conv) | numel(unit_conv) == 1, 'The number of elements in unit_conv must be 1 or equal to the number of cells in field.');
            if numel(unit_conv) == numel(field)
                for j = 1 : numel(field)
                    states{i}.(field{j}) = states{i}.(field{j})*unit_conv(j);
                end
            else
                for j = 1 : numel(field)
                    states{i}.(field{j}) = states{i}.(field{j})*unit_conv;
                end
            end
        end
        
    end
end

end

