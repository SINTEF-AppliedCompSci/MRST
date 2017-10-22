function [ states ] = changeUnits( states, fields, unit_conv )
%changeUnits converts the specified fields in states to the specified units. 
% 
% SYNOPSIS:
% 
%         [state] = changeUnits(state, fields, unit)
% 
% OUTPUT:
%  
%	state           -a fields with the values of state.fields changed to the specified units
% 
% REQUIRED PARAMETERS:
%
%	state           -The state structure produced by ChemicalModel/initstate. 
%                       State must be populated by the named variable before
%                       it can be retrived.
%
%	fields           -Name of the fields for which the unit change is to occur. 
%                       Can be a string or cell array of strings.
%
%	unit            -Numeric value of unit conversion. If numel(unit) == 1 
%                       and fields has multple entries then the one value of 
%                       unit will be applied to all fieldss. Otherwise if fields
%                       is a cell with multple entries then fields and unit must 
%                       have the same size, numel(unit) == numel(fields).unit_conv = unit_conv.^-1;
%
% EXAMPLE:
% 
%   change the units of a single fields
% 
%         [state] = changeUnits(state, 'activities', mol/litre)
% 
%   change the units of a multiple fieldss
% 
%         fieldss = {'acitivities', 'surfaceCharges'}
%         units = [mol/litre, mili*Coulumb/(nano*meter)^2];
% 
%         state = chem.getProps(state, fieldss, units);
% 
%   If state is a cell array of structures the function will loop over each cell. 
% 
%   The numeric value of the unit conversion must be defined. The is a bank 
%   of unit conversions within the release. 
% 

assert(isnumeric(unit_conv), 'The variable "unit_conv" must be numeric.');
unit_conv = unit_conv.^-1;

if numel(states) == 1
    
    if ~iscell(fields)
        states.(fields) = states.(fields)*unit_conv;
    else
        assert(numel(fields) == numel(unit_conv) | numel(unit_conv) == 1, 'The number of elements in unit_conv must be 1 or equal to the number of cells in fieldss.');
        if numel(unit_conv) == numel(fields)
            for i = 1 : numel(fields)
                states.(fields{i}) = states.(fields{i})*unit_conv(i);
            end
        else
            for i = 1 : numel(fields)
                states.(fields{i}) = states.(fields{i})*unit_conv;
            end
        end
    end
        
else

    for i = 1 : numel(states)
    
        if ~iscell(fields)
            states{i}.(fields) = states{i}.(fields)*unit_conv;
        else
            assert(numel(fields) == numel(unit_conv) | numel(unit_conv) == 1, 'The number of elements in unit_conv must be 1 or equal to the number of cells in fieldss.');
            if numel(unit_conv) == numel(fields)
                for j = 1 : numel(fields)
                    states{i}.(fields{j}) = states{i}.(fields{j})*unit_conv(j);
                end
            else
                for j = 1 : numel(fields)
                    states{i}.(fields{j}) = states{i}.(fields{j})*unit_conv;
                end
            end
        end
        
    end
end

end

