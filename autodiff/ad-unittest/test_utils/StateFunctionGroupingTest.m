function tests = StateFunctionGroupingTest
%Test suite for StateFunctionGrouping's dependency resolution and caching
%
% DESCRIPTION:
%   `StateFunctionGrouping`/`StateFunction` is the generic
%   dependency-graph property-evaluation framework underlying every
%   modern MRST model's discretization (flow/PVT property functions
%   etc.), but had no isolated unit test of its core behavior:
%   resolving a dependency before the function that needs it, and
%   caching evaluated values so repeated `get` calls do not
%   re-evaluate. This suite exercises that behavior directly with two
%   minimal mock state functions (`MockLeafStateFunction`, with no
%   dependencies, and `MockDependentStateFunction`, which depends on
%   the former), defined in `ad-unittest/setup/`.
%
% SEE ALSO:
%   StateFunctionGrouping, StateFunction, MockLeafStateFunction,
%   MockDependentStateFunction, functiontests

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

   tests = functiontests(localfunctions);
end

%--------------------------------------------------------------------------

function setupOnce(t)                                            %#ok<DEFNU>
   mrstModule reset
   mrstModule add ad-unittest ad-core
end

%--------------------------------------------------------------------------

function setup(t)                                                 %#ok<DEFNU>
   % Fresh model/grouping and reset evaluation counters for each test,
   % so counts from one test cannot leak into another.
   model = MockPhysicalModel();
   props = StateFunctionGrouping('props');
   props = props.setStateFunction('a', MockLeafStateFunction(model));
   props = props.setStateFunction('b', MockDependentStateFunction(model));

   t.TestData.model = model;
   t.TestData.props = props;
   t.TestData.state = props.initStateFunctionContainer(struct('x1', 5, 'x2', 3));

   MockLeafStateFunction.evaluationCount(-1);
   MockDependentStateFunction.evaluationCount(-1);
end

%--------------------------------------------------------------------------

function testDependencyIsResolvedBeforeDependent(t)              %#ok<DEFNU>
   % Requesting 'b' (which depends on 'a') must transparently evaluate
   % 'a' first and combine the values correctly: b = 2*x1 + x2.
   [bVal, state] = t.TestData.props.get(t.TestData.model, t.TestData.state, 'b'); %#ok<ASGLU>

   verifyEqual(t, bVal, 2*5 + 3);
   verifyEqual(t, MockLeafStateFunction.evaluationCount(0), 1);
   verifyEqual(t, MockDependentStateFunction.evaluationCount(0), 1);
end

%--------------------------------------------------------------------------

function testRepeatedGetIsCached(t)                               %#ok<DEFNU>
   % A second `get` of the same property on the same state object must
   % not trigger re-evaluation of either it or its dependency.
   [~, state] = t.TestData.props.get(t.TestData.model, t.TestData.state, 'b');
   [bVal2, ~] = t.TestData.props.get(t.TestData.model, state, 'b');

   verifyEqual(t, bVal2, 2*5 + 3);
   verifyEqual(t, MockLeafStateFunction.evaluationCount(0), 1);
   verifyEqual(t, MockDependentStateFunction.evaluationCount(0), 1);
end

%--------------------------------------------------------------------------

function testDependencyGetAfterDependentIsAlsoCached(t)          %#ok<DEFNU>
   % Once 'b' has triggered evaluation of its dependency 'a', a direct
   % `get` of 'a' afterwards must reuse the cached value rather than
   % evaluating it a second time.
   [~, state] = t.TestData.props.get(t.TestData.model, t.TestData.state, 'b');
   [aVal, ~]  = t.TestData.props.get(t.TestData.model, state, 'a');

   verifyEqual(t, aVal, 2*5);
   verifyEqual(t, MockLeafStateFunction.evaluationCount(0), 1);
end

%--------------------------------------------------------------------------

function testDirectDependencyGetDoesNotEvaluateDependent(t)      %#ok<DEFNU>
   % Requesting only the leaf property 'a' must not trigger evaluation
   % of the (unrelated, from 'a''s perspective) dependent property 'b'.
   [aVal, ~] = t.TestData.props.get(t.TestData.model, t.TestData.state, 'a');

   verifyEqual(t, aVal, 2*5);
   verifyEqual(t, MockLeafStateFunction.evaluationCount(0), 1);
   verifyEqual(t, MockDependentStateFunction.evaluationCount(0), 0);
end
