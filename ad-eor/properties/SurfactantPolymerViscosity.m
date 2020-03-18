classdef SurfactantPolymerViscosity < BlackOilViscosity
    methods
        function gp = SurfactantPolymerViscosity(model, varargin)
            gp@BlackOilViscosity(model, varargin{:});
            gp = addPropertyDependence(gp, {'PolymerViscMultiplier'});
            gp = addPropertyDependence(gp, {'SurfactantViscMultiplier'});
            assert(model.water);
        end
        
        function mu = evaluateOnDomain(prop, model, state)
            % The multipliers are only applied to the water phase
            mu   = prop.evaluateOnDomain@BlackOilViscosity(model, state);
            mPol = prop.getEvaluatedDependencies(state, 'PolymerViscMultiplier');
            mSft = prop.getEvaluatedDependencies(state, 'SurfactantViscMultiplier');
            ph   = model.getPhaseNames(); iW = find(ph=='W');
            mu{iW} = mu{iW} .* mPol .* mSft;
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
