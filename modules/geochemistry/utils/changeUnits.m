function [ states ] = changeUnits( states, fields, unit_conv )
% changeUnits converts the specified fields in states to the specified units. 
% 
% SYNOPSIS:
%   [state] = changeUnits(state, fields, unit)
%
% REQUIRED PARAMETERS:
%   state - The state structure produced by ChemicalModel/initstate. State 
%           must be populated by the named variable before they can be retrived.
%
%	fields - Name of the fields for which the unit change is to occur. Can 
%           be a string or cell array of strings.
%
%	unit - Numeric value of unit conversion. If numel(unit) == 1 and fields 
%           has multple entries then the one value of unit will be applied 
%           to all fields. Otherwise if fields is a cell with multple entries
%           then fields and unit must have the same size, 
%           numel(unit) == numel(fields).unit_conv = unit_conv.^-1;
% 
% RETURNS:
%   state - the state variable with the fields state.fields converted from
%           SI units to the units specified in units
%
% EXAMPLE:
%   % change the units of a single fields
% 
%         [state] = changeUnits(state, 'activities', mol/litre)
% 
%   %change the units of a multiple fields
% 
%         fields = {'acitivities', 'surfaceCharges'}
%         units = [mol/litre, mili*Coulumb/(nano*meter)^2];
% 
%         state = chem.getProps(state, fields, units);
% 
%   If state is a cell array of structures the function will loop over each cell. 
% 
%   The numeric value of the unit conversion must be defined. There is a bank 
%   of unit conversions within the release. 
% 

%{
Copyright 2009-2017 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
assert(isnumeric(unit_conv), 'The variable "unit_conv" must be numeric.');
unit_conv = unit_conv.^-1;

if numel(states) == 1
    
    if ~iscell(fields)
        states.(fields) = states.(fields)*unit_conv;
    else
        assert(numel(fields) == numel(unit_conv) | numel(unit_conv) == 1, 'The number of elements in unit_conv must be 1 or equal to the number of cells in fieldss.');
        if numel(unit_conv) == 1
            unit_conv = repmat(unit_conv, numel(fields), 1);
        end
        for i = 1 : numel(fields)
            f = fields{i};
            if iscell(states.(f))
                states.(f) = cellfun(@(x) x*unit_conv(i), states.(f), 'UniformOutput', false);
            else
                states.(f) = states.(f)*unit_conv(i);
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

