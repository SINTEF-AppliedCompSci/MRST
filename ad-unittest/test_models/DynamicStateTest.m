function tests = DynamicStateTest
    %Test suite for 'mergeOrderedArrays'
    %
    % SEE ALSO:
    %   functiontests, makeScheduleConsistent

    tests = functiontests(localfunctions);
end

%--------------------------------------------------------------------------

function setupOnce(t)
    mrstModule add ad-props ad-blackoil ad-core
    G = cartGrid(10); G = computeGeometry(G);
    rock = makeRock(G, 1, 1);
    fluid = initSimpleADIFluid();
    model = ReservoirModel(G, rock, fluid);
    t.TestData.model = model;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Constructor tests                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Retrieval tests                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function testDynamicRetrieval(t)
    state = struct();
    
    variables = {1, 2, 3};
    names = {'a', 'b', 'c'};
    
    dyn = DynamicState(state, names, variables);
    b = dyn.a;
    
    verifyEqual(t, variables{1}, b)
end

function testDynamicRetrievalInDynVars(t)
    state = struct('d', (1:10)');
    variables = {(5:10)', (2:10)', (3:10)'};
    names = {'a', 'b', 'c'};
    
    dyn = DynamicState(state, names, variables);
    b = dyn.a(:, 1);
    
    verifyEqual(t, variables{1}, b)
    
    subs = 5:6;
    b = dyn.b(subs);
    verifyEqual(t,  variables{2}(subs), b)
    
    b = dyn.b(subs, :);
    verifyEqual(t,  variables{2}(subs, :), b)
end


function testDynamicRetrievalInStaticVars(t)
    state = struct('d', (1:10)');
    variables = {(5:10)', (2:10)', (3:10)'};
    names = {'a', 'b', 'c'};
    
    dyn = DynamicState(state, names, variables);
    b = dyn.d(:, 1);
    
    verifyEqual(t, state.d, b)
    
    subs = 5:6;
    b = dyn.d(subs);
    verifyEqual(t, state.d(subs), b)
    
    b = dyn.d(subs, :);
    verifyEqual(t, state.d(subs, :), b)
end

% Test scalar variable
function testGetPropPressure(t)
    model = t.TestData.model;
    state = initResSol(model.G, 100);
    
    dyn = DynamicState(state);
    result = model.getProp(dyn, 'pressure');
    verifyEqual(t, state.pressure, result)
end

function testGetPropPressureAD(t)
    model = t.TestData.model;
    state = initResSol(model.G, 100);
    
    p = model.AutoDiffBackend.initVariablesAD(state.pressure);
    
    dyn = DynamicState(state, {'pressure'}, {p});
    result = model.getProp(dyn, 'pressure');
    verifyEqual(t, p, result)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Setting variables                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
