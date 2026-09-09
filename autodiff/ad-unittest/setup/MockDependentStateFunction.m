classdef MockDependentStateFunction < StateFunction
    % Minimal single-dependency StateFunction (depends on the state
    % function named 'a' in the same grouping) used to test
    % StateFunctionGrouping's dependency resolution and caching in
    % StateFunctionGroupingTest. Counts how many times it has actually
    % been evaluated, so tests can assert on cache hits/misses.
    methods
        function prop = MockDependentStateFunction(model)
            prop = prop@StateFunction(model);
            prop = prop.dependsOn('a');
        end

        function v = evaluateOnDomain(prop, model, state) %#ok<INUSL>
            MockDependentStateFunction.evaluationCount(1);
            a = prop.getEvaluatedDependencies(state, 'a');
            v = a + state.x2;
        end
    end

    methods (Static)
        function n = evaluationCount(bump)
            persistent count
            if isempty(count)
                count = 0;
            end
            if nargin > 0
                if bump < 0
                    count = 0; % reset
                else
                    count = count + bump;
                end
            end
            n = count;
        end
    end
end

%{
Copyright 2009-2026 SINTEF Digital, Mathematics & Cybernetics.

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
