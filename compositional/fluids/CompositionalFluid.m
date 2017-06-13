classdef CompositionalFluid
    properties
        % Critical temperature in Kelvin
        Tcrit
        % Critical pressure in Pascal
        Pcrit
        % Critical volume in m^3 / mol
        Vcrit
        % Acentric factors (dimensionless)
        acentricFactors
        % Component mass (kg / mol)
        molarMass
        T_ref
        P_ref
        % Names of each component. Must be unique.
        names
    end
    
    properties ( Access = private)
        % Binary interaction coefficients
        bic
    end
    
    methods
        function fluid = CompositionalFluid(names, Tcrit, Pcrit, Vcrit, acentricFactors, molarMass, T_ref, P_ref)
            if nargin == 0
                disp('Empty fluid - no validation');
                return
            end
            
            fluid.names = names;
            fluid.Tcrit = Tcrit;
            fluid.Pcrit = Pcrit;
            fluid.Vcrit = Vcrit;
            fluid.acentricFactors = acentricFactors;
            fluid.molarMass = molarMass;
            fluid.T_ref = T_ref;
            fluid.P_ref = P_ref;
            ncomp = numel(fluid.names);
            
            fluid.bic = zeros(ncomp, ncomp);
        end
        
        function n = getNumberOfComponents(fluid)
            n = numel(fluid.names);
        end
        
        function bic = getBinaryInteraction(fluid)
            bic = fluid.bic;
        end
        
        function fluid = setBinaryInteraction(fluid, input)
            % Set BIC via a matrix. Must be symmetric and ncomp by ncomp
            ncomp = fluid.getNumberOfComponents();
            assert(size(input, 1) == ncomp)
            assert(size(input, 1) == size(input, 2));
            assert(all(all(input == input')));
            fluid.bic = input;
        end
    end
    
end

%{
Copyright 2009-2017 SINTEF ICT, Applied Mathematics.

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
