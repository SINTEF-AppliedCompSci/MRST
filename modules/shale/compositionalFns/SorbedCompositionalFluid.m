classdef SorbedCompositionalFluid < TableCompositionalMixture
    % Create MRST compositional from stored tables (taken from CoolProp)
    properties
        isotherm % isotherm parameters describing fluid in the nanopores
    end
    
    methods
        function fluid = SorbedCompositionalFluid(names,isotherm)
            % Create fluid from names
            %
            % SYNOPSIS:
            %   f = NanoCompositionalFluid({'Methane', 'n-Decane'});
            %
            % PARAMETERS:
            %   names - Cell array of valid names. See `getFluidList` for
            %           valid names.
            % RETURNS:
            %   fluid - Initialized fluid.
            fluid = fluid@TableCompositionalMixture(names);
            %OMO Edit: incorporate the isotherm into the fluid object
            fluid.isotherm = isotherm;
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
