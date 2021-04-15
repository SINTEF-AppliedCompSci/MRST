classdef CoolPropsCompositionalMixture < CompositionalMixture
    % Create MRST compositional fluid by using the open CoolProp library
    % (can be downloaded from http://www.coolprop.org/)
    properties
        
    end
    
    methods
        function fluid = CoolPropsCompositionalMixture(names)
            % Create fluid from names
            %
            % SYNOPSIS:
            %   f = TableCompositionalMixture({'Methane', 'n-Decane'});
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
            
            validChoices = CoolPropsCompositionalMixture.getFluidList();
            ok = ismember(lower(names), lower(validChoices));
            
            if ~all(ok)
                s ='Unable to create fluid. The following names were not known to coolprops: ';
                msg = [s, sprintf('%s ', names{~ok})];
                error(msg);
            end
            
            for i = 1:numel(names)
                isF = strcmpi(names{i}, validChoices);
                n = validChoices{isF};
                Tcrit(i) = CoolProp.Props1SI(n, 'Tcrit');
                Pcrit(i) = CoolProp.Props1SI(n, 'Pcrit');
                rhocrit(i) = CoolProp.Props1SI(n, 'rhocrit');
                acc(i) = CoolProp.Props1SI(n, 'acentric');
                molarMass(i) = CoolProp.Props1SI(n, 'molarmass');
            end
            assert(all(isfinite(Tcrit)));
            assert(all(isfinite(Pcrit)));
            assert(all(isfinite(rhocrit)));
            assert(all(isfinite(acc)));
            assert(all(isfinite(molarMass)));
            Vcrit = molarMass./rhocrit;
            fluid = fluid@CompositionalMixture(names, Tcrit, Pcrit, Vcrit, acc, molarMass);
        end
    end
    methods (Static)
        function varargout = getFluidList()
            names = num2str(CoolProp.get_global_param_string('FluidsList'));
            names = strsplit(names, ',');
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
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

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
