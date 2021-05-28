classdef SurfactantAdsorption < StateFunction
% Surfactant adsorption. The quantity of adsorped surfactant will affect the
% mass conservation equation for surfactant.
    properties
    end

    methods
        function gp = SurfactantAdsorption(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'surfactant', 'surfactantmax'}, 'state'); % check mechanism
            assert(model.water);
            assert(isfield(model.fluid,'adsInxSft'));
        end

        function ads = evaluateOnDomain(prop, model, state)
            [cs, csmax] = model.getProps(state, 'surfactant', 'surfactantmax');
            if model.fluid.adsInxSft == 2
                ce = max(cs, csmax);
            else
                ce = cs;
            end
            ads = prop.evaluateFluid(model, 'surfads', ce);
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
