classdef PhaseMultipliers < StateFunction
% This class is a container class which is used to store and compute the product
% of multipliers that are meant to be apply to a given quantity. An example is
% the ViscosityMultipliers container which is set in the
% ThreePhaseSurfactantPolymerModel.
%
% Each multiplier is defined as a statefunction, see PolymerEffViscMult for
% an example of multiplier.
    
    properties
        multNames; % contains the name of the multipliers.
        operator;  % the default operator is multiplication but an other
                   % operator can be given here
    end
    
    methods
        function M = PhaseMultipliers(model, varargin)
            M@StateFunction(model, varargin{:});
            M.label = 'M';
        end

        function prop = addMultiplier(prop, model, name, phase)
            phind = model.getPhaseIndex(phase);
            nph = model.getNumberOfPhases();
            multName = cell(1, nph);
            multName{phind} = name;
            multNames = prop.multNames;
            multNames{end + 1} = multName;
            prop.multNames = multNames;
            prop = prop.dependsOn(name);
        end
        
        
        function M = evaluateOnDomain(prop, model, state)
            nph = model.getNumberOfPhases();
            multNames = prop.multNames;
            nmult = numel(multNames);
            
            M = cell(1, nph);
            for ph = 1:nph
                mult = 1;
                % activemult will be set to true if there exists an active multiplier for this phase
                activemult = false;
                for i = 1 : nmult
                    multName = multNames{i};
                    if ~isempty(multName{ph})
                        activemult = true;
                        multi = prop.getEvaluatedDependencies(state, multName{ph});
                    if isempty(prop.operator)
                            mult = mult .* multi;
                    else
                            mult = prop.operator(mult, multi);
                        end
                    end
                end
                if activemult
                    M{ph} = mult;
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
