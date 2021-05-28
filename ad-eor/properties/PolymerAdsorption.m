classdef PolymerAdsorption < StateFunction
% Polymer adsorption. The quantity of adsorped polymer will affect the mass
% conservation equation for polymer and the relative permeability of water,
% see PolymerPermReduction.
    
    properties
    end

    methods
        function gp = PolymerAdsorption(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'polymer', 'polymermax'}, 'state'); % check mechanism
            assert(model.water && isfield(model.fluid, 'adsInx'));
        end

        function ads = evaluateOnDomain(prop, model, state)
            [cp, cpmax] = model.getProps(state, 'polymer', 'polymermax');
            if model.fluid.adsInx == 2
                ce = max(cp, cpmax);
            else
                ce = cp;
            end
            ads = prop.evaluateFluid(model, 'ads', ce);
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
