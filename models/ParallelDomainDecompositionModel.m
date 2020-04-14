classdef ParallelDomainDecompositionModel < DomainDecompositionModel
   
    properties
        localModel
        pool = [];
        domains = [];
        NumWorkers;
    end
    
    methods
        
        %-----------------------------------------------------------------%
        function [state, report, extra] = solveSubDomains(model, state0, dt, drivingForces, state)
            
            % Initialize
            [initState, startState] = deal(state);
            initState = stripState(initState);
            state0    = stripState(state0);
            sds = model.subdomainSetup;
            lm  = model.localModel;
            mapStateOnWorkers = true;
            if mapStateOnWorkers
                add = @(in1, in2) addStates(in1, in2, initState);
            end
            ticBytes(model.pool);
            addProblems = @(p1, p2) model.addProblems(p1, p2);
            spmd
                localTimer = tic(); % Time the local solve (db only)
                map = [];
                [pFull, pSub] = deal([]);
                finalState = initState;
                iterations = zeros(lm.G.cells.num,1);
                for i = 1:numel(lm.domains)
                    % Get cell subset
                    cells = lm.partition == lm.domains(i);
                    % Solve subdomain
                    [state, substate, subreport, mappings, subproblems] = lm.solveSubDomain(sds{i}, state0, dt, drivingForces, initState);
                    if lm.outputJacobians
                        % Compute Jacobians for subproblem
                        pFull = lm.addProblems(pFull, subproblems.full);
                        pSub  = lm.addProblems(pSub , subproblems.sub);
                    end
                    % Update iterations
                    iterations(cells) = subreport.Iterations;
                    finalState = mapState(finalState, substate, mappings);
                    map = mergeMappings(map, mappings);
                end
                t_localSolve = toc(localTimer);
                % Sum Jacobians and gather on worker 1
                pFull = gop(addProblems, pFull, 1);
                pSub  = gop(addProblems, pSub, 1);
                % Sum iterations and gather on worker 1
                iterations = gop(@plus, iterations, 1);
                if numel(lm.domains) > 0
                    finalState = getSubState(finalState, map);
                end
                if mapStateOnWorkers
                    in = struct('mappings', map, 'state', finalState);
                    out = gop(add, in, 1);
                    if labindex == 1
                        finalState = out.state;
                    end
                end
                t_localComm = toc(localTimer) - t_localSolve; %#ok
            end
            bytes = tocBytes(model.pool); %#ok
            iterations = iterations{1};
            timer = tic();
            if mapStateOnWorkers
                state = finalState{1};
                state.iterations = startState.iterations;
            else
                state = startState;
                for i = 1:model.NumWorkers
                    if ~isempty(map{i})
                        state = mapState(state, finalState{i}, map{i});
                    end
                end
            end
            if model.outputJacobians
                problems = struct('full', pFull{1}, 'sub', pSub{1});
                extra     = struct('subproblems', problems, 'initState', startState);
            else
                extra = [];
            end
            t_comm = toc(timer); %#ok
            % Make step report
            report = model.makeDomainStepReport('Iterations', iterations);
            
        end
        
        %-----------------------------------------------------------------%
        function model = validateModel(model, varargin)
            if isempty(model.pool)
                model.pool = gcp();
                model.NumWorkers = model.pool.NumWorkers;
            end
            model = validateModel@WrapperModel(model, varargin{:});
            m = model;
            m.G = [];
            rmodel = m.getReservoirModel();
            if isfield(rmodel, 'FacilityModel')
                rmodel.FacilityModel = [];
            end
            spmd
               lm = m;
               lm.G = lm.parentModel.G;
%                lm.parentModel.FacilityModel.ReservoirModel = lm.parentModel;
            end
            model.localModel = lm;
            model = model.updateSubdomainSetup(model.partition);
        end
        
        %-----------------------------------------------------------------%
        function model = updateSubdomainSetup(model, partition)
            lm = model.localModel;
            spmd
                lm.partition = partition;
                dom = getLocalDomains(lm.NumWorkers, lm.partition, labindex);
                setup = cell(1,numel(dom));
                for i = 1:numel(dom)
                    cells = partition == dom(i);
                    setup{i} = lm.getSubdomainSetup(cells);
                    G = setup{i}.Model.G;
                    G = rmfield(G, 'nodes');
                    setup{i}.Model.G = G;
                    setup{i}.Model.parentModel.G = G;
                end
                lm.domains = dom;
            end
            model.localModel = lm;
            model.subdomainSetup = setup;
        end
        
    end
    
end
%-------------------------------------------------------------------------%
function state = stripState(state)
    fields = {'sMax', 'iterations', 'FacilityState', 'K', ...
        'eos', 'FacilityState', 'switched', 'switchCount', 'dpRel', 'dpAbs'};
%     fields = {'flux', 'sMax', 'iterations', 'FacilityState'};
    for f = fields
        if isfield(state, f{1})
            state = rmfield(state, f{1});
        end
    end
end

%-------------------------------------------------------------------------%
function domains = getLocalDomains(numWorkers, p, i)
    nDomains = repmat(ceil(max(p)/numWorkers), 1, numWorkers);
    extra    = sum(nDomains) - max(p);
    nDomains(end-extra+1:end) = nDomains(end-extra+1:end) - 1;
    domains = [0, cumsum(nDomains)] + 1;
    domains = domains(i):domains(i+1)-1;
end

%-------------------------------------------------------------------------%
function mappings = mergeMappings(mappings1, mappings2)
    if isempty(mappings1)
        mappings1 = mappings2;
    end
    if islogical(mappings1) && islogical(mappings2)
        mappings = mappings1 | mappings2;
        return
    elseif iscell(mappings1) && iscell(mappings2)
        mappings = union(mappings1, mappings2);
        return
    elseif isnumeric(mappings1) && isnumeric(mappings2)
        mappings = mappings1;
        return
    end
    flds = fieldnames(mappings1);
    for i = 1:numel(flds)
        map1 = mappings1.(flds{i});
        map2 = mappings2.(flds{i});
        map  = mergeMappings(map1, map2);
        mappings.(flds{i}) = map;
    end
end

%-------------------------------------------------------------------------%
function out = addStates(in1, in2, dummyState)

    map1 = in1.mappings;
    map2 = in2.mappings;
    if isempty(map1) || isempty(map2)
        if isempty(map2)
            out = in1;
        else
            out = in2;
        end
        return
    end
    substate1 = in1.state;
    substate2 = in2.state;
    dummyState = mapState(dummyState, substate1, map1);
    dummyState = mapState(dummyState, substate2, map2);
    map = mergeMappings(map1, map2);
    substate = getSubState(dummyState, map);
    out = struct('state', substate, 'mappings', map);
    
end