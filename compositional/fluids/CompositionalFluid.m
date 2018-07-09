classdef CompositionalFluid
    properties
        Tcrit % Critical temperature in Kelvin
        Pcrit % Critical pressure in Pascal
        Vcrit % Critical volume in m^3 / mol
        acentricFactors % Acentric factors (dimensionless)
        molarMass % Component mass (kg / mol)
        names % Names of each component. Each name must be unique.
    end
    
    properties ( Access = protected)
        bic % Binary interaction coefficients
    end
    
    methods
        function fluid = CompositionalFluid(names, Tcrit, Pcrit, Vcrit, acentricFactors, molarMass, varargin)
            if nargin == 0
                disp('Empty fluid - no validation');
                return
            end
            
            for i = 1:numel(names)
                cts = sum(strcmp(names{i}, names));
                assert(cts == 1, ...
                    ['Component ', num2str(i), ': ', names{i}, ' occurs multiple times. Components must be unique.']);
            end
            fluid.names = names;
            fluid.Tcrit = Tcrit;
            fluid.Pcrit = Pcrit;
            fluid.Vcrit = Vcrit;
            fluid.acentricFactors = acentricFactors;
            fluid.molarMass = molarMass;
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
            
            if isvector(input)
                n_el = ncomp*(ncomp-1)/2;
                assert(numel(input) == n_el);
                
                tmp = zeros(ncomp, ncomp);
                offset = 0;
                for i = 2:ncomp
                    cts = i-1;
                    tmp(i, 1:i-1) = input(offset + (1:cts));
                    offset = offset + cts;
                end
                input = tmp + tmp';
            end
            assert(size(input, 1) == ncomp)
            assert(size(input, 1) == size(input, 2));
            assert(all(all(input == input')));
            fluid.bic = input;
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
