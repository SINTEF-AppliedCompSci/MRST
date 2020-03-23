classdef MultiplierContainer < StateFunction
% This class is a container class which is used to store and compute the product
% of multipliers that are meant to be apply to a given quantity. An example is
% the ViscosityMultipliers container which is set in the
% ThreePhaseSurfactantPolymerModel.
%
% Each multiplier is defined as a statefunction, see PolymerViscMultiplier for
% an example of multiplier.
    
    properties
        multNames; % contains the name of the multipliers.
        operator;  % the default operator is multiplication but an other
                   % operator can be given here
    end
    
    methods
        function M = MultiplierContainer(model, varargin)
            M@StateFunction(model, varargin{:});
            M.label = 'M';
        end

        function prop = addMultiplier(prop, model, name, phase)
            phind = model.getPhaseIndex(phase);
            nph = model.getNumberOfPhases();
            multname = cell(1, nph);
            multName{phind} = name;
            multNames = prop.multNames;
            multNames{end + 1} = multName;
            prop.multNames = multNames;
            prop = prop.dependsOn(name);
        end
        
        
        function M = evaluateOnDomain(prop, model, state)
            nph = numel(prop.multNames);
            M = cell(1, nph);
            for ph = 1:nph
                multPhaseNames = prop.multNames{ph};
                if numel(multPhaseNames) == 0
                    continue
                end
                mult = 1;
                for i = 1:numel(multPhaseNames)
                    % Multiplier can be compund, i.e. M = A*B*C
                    mi = prop.getEvaluatedDependencies(state, multPhaseNames{i});
                    if isempty(prop.operator)
                        mult = mult .* mi;
                    else
                        mult = prop.operator(mult, mi);
                    end
                end
                M{ph} = mult;
            end
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
