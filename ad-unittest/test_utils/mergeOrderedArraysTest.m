function tests = mergeOrderedArraysTest
    %Test suite for 'mergeOrderedArrays'
    %
    % SEE ALSO:
    %   functiontests, makeScheduleConsistent

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

function setupOnce(t)

end

%--------------------------------------------------------------------------
function testVaried(t)
    a = [1; 3; 4];
    b = [1; 2; 3; 4];
    c = mergeOrderedArrays(a, b);

    answer = [1; 2; 3; 4];
    verifyEqual(t, c, answer)
end


function testOrderedEqual(t)
    a = [1; 2; 3];
    b = [1; 2; 3];
    c = mergeOrderedArrays(a, b);

    verifyEqual(t, a, c)
end


function testLastElement(t)
    a = [1; 2];
    b = [1; 2; 3];
    c = mergeOrderedArrays(a, b);

    answer = [1; 2; 3];
    verifyEqual(t, c, answer)
end
