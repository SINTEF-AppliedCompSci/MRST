classdef SurfactantPolymerViscosity < BlackOilViscosity
    methods
        function gp = SurfactantPolymerViscosity(prop, varargin)
            gp@BlackOilViscosity(prop, varargin{:});
            gp = addPropertyDependence(gp, {'PolymerViscMultiplier'});
            gp = gp.dependsOn('surfactant', 'state');
            gp = gp.dependsOn('pressure', 'state');
        end
        
        function mu = evaluateOnDomain(prop, model, state)
            fluid = model.fluid;
            p = model.getProps(state, 'pressure');
            cs = model.getProps(state, 'surfactant');
            muWMultp = prop.getEvaluatedDependencies(state, 'PolymerViscMultiplier');
            muWMults = fluid.muW(p)/fluid.muWr;
            mu = prop.evaluateOnDomain@BlackOilViscosity(model, state);
            mu{1} = model.fluid.muWSft(cs);
            mu{1} = mu{1}.*muWMults.*muWMultp;            
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
