classdef PhaseMultipliers < StateFunction
    properties
        multNames
    end
    
    methods
        function M = PhaseMultipliers(model, names, varargin)
            M@StateFunction(model, varargin{:});
            nph = model.getNumberOfPhases();
            assert(numel(names) == nph);
            M.multNames = cell(1, nph);
            for ph = 1:nph
                n = names{ph};
                if isempty(n)
                    n = {};
                elseif ischar(n)
                    n = {n};
                end
                assert(iscell(n));
                M = M.dependsOn(n);
                M.multNames{ph} = n;
            end
            M.label = 'M';
        end
        
        function M = evaluateOnDomain(prop, model, state)
            nph = numel(prop.multNames);
            M = cell(1, nph);
            for ph = 1:nph
                phaseNames = prop.multNames{ph};
                if numel(phaseNames) == 0
                    continue
                end
                mult = 1;
                for i = 1:numel(phaseNames)
                    % Multiplier can be compund, i.e. M = A*B*C
                    mi = prop.getEvaluatedDependencies(state, phaseNames{i});
                    mult = mult * mi;
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
