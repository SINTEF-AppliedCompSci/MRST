classdef BlackOilSurfactantViscosity < BlackOilViscosity
    properties
    end

    methods
        function gp = BlackOilSurfactantViscosity(model, varargin)
            gp@BlackOilViscosity(model, varargin{:});
            gp = gp.dependsOn('PhasePressures');
            assert(model.water);
            assert(all(isfield(model.fluid,{'muW','muWR'})));
        end

        function mu = evaluateOnDomain(prop, model, state)
           cs = model.getProps(state, 'surfactant');
           p = model.getProps(state, 'PhasePressures');
           pW = p{1};
           mu = prop.evaluateOnDomain@BlackOilViscosity(model, state);
           mu{1} = model.fluid.muWSft(cs);
           muWMultPressure = model.fluid.muW(pW)/model.fluid.muWr;
           
           mu{1} = mu{1}.*muWMultPressure;
        end
    end
end

%{
Copyright 2009-2019 SINTEF Digital, Mathematics & Cybernetics.

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
