classdef SurfactantViscMultiplier < StateFunction
% Surfactant viscosity multiplier. This value is multiplied to the viscosity and
% account in this way for the effect of surfactant on the viscosity.
    
    properties
    end

    methods
        function gp = SurfactantViscMultiplier(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'surfactant'}, 'state'); % check mechanism
            assert(model.water && all(isfield(model.fluid,{'muWr','muWSft'})));
        end

        function mSft = evaluateOnDomain(prop, model, state)
            % The multiplier is only applied to the water phase
            cs   = model.getProps(state,'surfactant');
            mSft = prop.evaluateFluid(model, 'muWSft', cs);
            mSft = mSft / model.fluid.muWr;
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
