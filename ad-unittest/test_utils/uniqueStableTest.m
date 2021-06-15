function tests = uniqueStableTest
%Test suite for 'uniqueStable' wrapper
%
% SEE ALSO:
%   uniqueStable, functiontests

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

   tests = functiontests(localfunctions);
end

%--------------------------------------------------------------------------

function setupOnce(t)                                           %#ok<DEFNU>
   t.TestData.a = [ 3, 3, 3 ; ...
                    1, 4, 1 ; ...
                    2, 2, 2 ; ...
                    1, 1, 1 ; ...
                    2, 2, 2 ; ...
                    3, 3, 3 ];
end

%--------------------------------------------------------------------------

function testStableElemSort(t)                                  %#ok<DEFNU>
   c = uniqueStable(t.TestData.a, 'use_fallback');
   e = [ 3 ; 1 ; 2 ; 4 ];  % Expected solution

   verifyEqual(t, c, e)
end

%--------------------------------------------------------------------------

function testStableElemIA(t)                                    %#ok<DEFNU>
   [c, ia] = uniqueStable(t.TestData.a, 'use_fallback');

   verifyEqual(t, c, t.TestData.a(ia));
end

%--------------------------------------------------------------------------

function testStableElemIC(t)                                    %#ok<DEFNU>
   [c, ic, ic] = uniqueStable(t.TestData.a, 'use_fallback');    %#ok<ASGLU>

   verifyEqual(t, reshape(c(ic), size(t.TestData.a)), t.TestData.a);
end

%--------------------------------------------------------------------------

function testStableRowSort(t)                                   %#ok<DEFNU>
   c  = uniqueStable(t.TestData.a, 'rows', 'use_fallback');
   e1 = [ 3 ; 1 ; 2 ; 1 ];
   et = t.TestData.a([1, 2, 3, 4], :);

   verifyEqual(t, c(2,2), 4 );
   verifyEqual(t, c(:,1), e1);
   verifyEqual(t, c     , et);
end

%--------------------------------------------------------------------------

function testStableRowIA(t)                                     %#ok<DEFNU>
   [c, ia] = uniqueStable(t.TestData.a, 'rows', 'use_fallback');

   verifyEqual(t, c, t.TestData.a(ia, :));
end

%--------------------------------------------------------------------------

function testStableRowIC(t)                                     %#ok<DEFNU>
   [c, ic, ic] = uniqueStable(t.TestData.a, 'rows', ...
                              'use_fallback');                  %#ok<ASGLU>

   verifyEqual(t, c(ic,:), t.TestData.a);
end

%--------------------------------------------------------------------------

function testStableRowsOnColumn(t)                              %#ok<DEFNU>
   % Check that 'rows' on column vector is equivalent to no 'rows' at all

   e = [ 3 ; 4 ; 2 ; 1 ];

   [c, ia, ic] = uniqueStable(t.TestData.a(:,2), 'rows', 'use_fallback');
   verifyEqual(t, c    , e);
   verifyEqual(t, c    , t.TestData.a(ia,2));
   verifyEqual(t, c(ic), t.TestData.a(: ,2));

   [c2, ia2, ic2] = uniqueStable(t.TestData.a(:,2),      'use_fallback');
   verifyEqual(t, c , c2 );
   verifyEqual(t, ia, ia2);
   verifyEqual(t, ic, ic2);
end
