function tests = simpleEquilibriumTest
%Test suite for 'simpleEquilibrium'
%
% SEE ALSO:
%   simpleEquilibrium, functiontests

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

function setupOnce(t)                                           %#ok<DEFNU>
   G = cartGrid([10, 1, 10], [1, 1, 1]);
   t.TestData.G = computeGeometry(G);
   t.TestData.contacts = [.37, .8];
end

%--------------------------------------------------------------------------

function testOutputShape(t)                                     %#ok<DEFNU>
   G = t.TestData.G;
   s = simpleEquilibrium(G, t.TestData.contacts);

   verifySize(t, s, [G.cells.num, numel(t.TestData.contacts) + 1]);
end

%--------------------------------------------------------------------------

function testSaturationsInUnitRangeAndSumToOne(t)                %#ok<DEFNU>
   G = t.TestData.G;
   s = simpleEquilibrium(G, t.TestData.contacts);

   verifyTrue(t, all(s(:) >= -sqrt(eps) & s(:) <= 1 + sqrt(eps)));
   verifyEqual(t, sum(s, 2), ones(G.cells.num, 1), 'AbsTol', sqrt(eps));
end

%--------------------------------------------------------------------------

function testNoContactsGivesFullSaturation(t)                    %#ok<DEFNU>
   % With zero contacts, every cell is fully within the single (only)
   % phase, so the sole output column must be all ones.
   G = t.TestData.G;
   s = simpleEquilibrium(G, []);

   verifySize(t, s, [G.cells.num, 1]);
   verifyEqual(t, s, ones(G.cells.num, 1), 'AbsTol', sqrt(eps));
end

%--------------------------------------------------------------------------

function testAlternateDirectionVector(t)                         %#ok<DEFNU>
   % Passing an explicit direction vector must still produce a valid,
   % normalized saturation field of the same shape as the default case.
   G = t.TestData.G;
   contacts = t.TestData.contacts;

   s = simpleEquilibrium(G, contacts, [1, 0, 1]);

   verifySize(t, s, [G.cells.num, numel(contacts) + 1]);
   verifyEqual(t, sum(s, 2), ones(G.cells.num, 1), 'AbsTol', sqrt(eps));
end
