function tests = DynamicStateTest
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
function testConstructorOnlyState(t)
    state = struct();
    
    dyn = DynamicState(state);
end

function testConstructorVariables(t)
    state = struct();
    
    variables = {1, 2, 3};
    names = {'a', 'b', 'c'};
    
    dyn = DynamicState(state, variables, names);
end

function testDynamicRetrieval(t)
    state = struct();
    
    variables = {1, 2, 3};
    names = {'a', 'b', 'c'};
    
    dyn = DynamicState(state, names, variables);
    b = dyn.a;
    
    verifyEqual(t, variables{1}, b)
end

function testDynamicRetrievalInDynVars(t)
    state = struct('d', 5);
    
    variables = {1, 2, 3};
    names = {'a', 'b', 'c'};
    
    dyn = DynamicState(state, names, variables);
    b = dyn.a(:, 1);
    
    verifyEqual(t, variables{1}, b)
end


function testDynamicRetrievalInStaticVars(t)
    state = struct('d', 5);
    
    variables = {1, 2, 3};
    names = {'a', 'b', 'c'};
    
    dyn = DynamicState(state, names, variables);
    b = dyn.d(:, 1);
    
    verifyEqual(t, state.d, b)
end