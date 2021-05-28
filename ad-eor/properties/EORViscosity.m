classdef EORViscosity < Viscosity
% Viscosity for EOR model which is based on multipliers
% 
% This viscosity function is defined by multiplying the black-oil viscosity with
% multipliers where each multiplier correspond to a given EOR process. See
% SurfactantViscMultiplier and PolymerEffViscMult.
%   
% The statefunction ViscosityMultipliers collects all the multipliers and is
% setup by the model as a container (see helper class PhaseMultipliers).
% 
    
    methods
        function mu = EORViscosity(model, varargin)
            mu@Viscosity(model, varargin{:});
            mu = mu.dependsOn({'ViscosityMultipliers', 'BaseViscosity'});
        end
        
        function mu = evaluateOnDomain(prop, model, state)
            mu = prop.getEvaluatedDependencies(state, 'BaseViscosity');
            mult = prop.getEvaluatedDependencies(state, 'ViscosityMultipliers');
            for i = 1 : numel(mult)
                m = mult{i};
                if ~isempty(m)
                    mu{i} = mu{i}.*m;
                end
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
