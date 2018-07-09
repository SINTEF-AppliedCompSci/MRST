classdef TableCompositionalFluid < CompositionalFluid
    % Create MRST compositional from stored tables (taken from CoolProp)
    properties
        
    end
    
    methods
        function fluid = TableCompositionalFluid(names)
            % Create fluid from names
            %
            % SYNOPSIS:
            %   f = TableCompositionalFluid({'Methane', 'n-Decane'});
            %
            % PARAMETERS:
            %   names - Cell array of valid names. See `getFluidList` for
            %           valid names.
            % RETURNS:
            %   fluid - Initialized fluid.
            if ischar(names)
                names = {names};
            end
            
            ncomp = numel(names);
            [Tcrit, Pcrit, rhocrit, acc, molarMass] = deal(zeros(1, ncomp));
            
            fluids = coolPropFluidsStructs();
            validChoices = TableCompositionalFluid.getFluidList();
            
            ok = ismember(lower(names), lower(validChoices));
            
            if ~all(ok)
                s ='Unable to create fluid. The following names were not known to CoolProps table: ';
                msg = [s, sprintf('%s ', names{~ok})];
                error(msg);
            end
            
            for i = 1:numel(names)
                isF = strcmpi(names{i}, validChoices);
                str = fluids(isF);
                Tcrit(i) = str.Tcrit;
                Pcrit(i) = str.Pcrit;
                rhocrit(i) = str.rhocrit;
                acc(i) = str.acentric;
                molarMass(i) = str.molarmass;
            end
            Vcrit = molarMass./rhocrit;
            fluid = fluid@CompositionalFluid(names, Tcrit, Pcrit, Vcrit, acc, molarMass);
        end
    end
    methods (Static)
        function varargout = getFluidList()
            fluids = coolPropFluidsStructs();
            names = {fluids.name};
            if nargout > 0
                varargout{1} = names;
            else
                disp('Possible fluid choices are:');
                fprintf('%s\n', names{:});
            end
        end 
    end
end

%{
Copyright 2009-2018 SINTEF ICT, Applied Mathematics.

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
