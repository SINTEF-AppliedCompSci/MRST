classdef  EORRelativePermeability < BaseRelativePermeability
% Relative permeability for EOR models which are based on a multiplier approach
% and a generic black-oil relative permeability model
% 
% This relative function is defined by multiplying the black-oil rel.perm model
% with multipliers where each multiplier correspond to a given EOR process. See
% PolymerPermReduction.
%   
% The statefunction RelativePermeabilityMultipliers collects all the multipliers
% and is setup by the model as a container (see helper class
% PhaseMultipliers).
% 

    properties

    end

    methods
        function gp = EORRelativePermeability(model, varargin)
            gp@BaseRelativePermeability(model, varargin{:});
            gp = gp.dependsOn({'RelativePermeabilityMultipliers', 'BaseRelativePermeability'});
        end

        function kr = evaluateOnDomain(prop, model, state)
            kr   = prop.getEvaluatedDependencies(state, 'BaseRelativePermeability');
            mult = prop.getEvaluatedDependencies(state, 'RelativePermeabilityMultipliers');
            for i = 1 : numel(mult)
                m = mult{i};
                if ~isempty(m)
                    kr{i} = kr{i}.*m;
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
