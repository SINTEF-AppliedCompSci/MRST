function tests = mergeOrderedArraysTest
    %Test suite for 'mergeOrderedArrays'
    %
    % SEE ALSO:
    %   functiontests, makeScheduleConsistent

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
