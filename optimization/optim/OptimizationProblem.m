classdef OptimizationProblem < BaseEnsemble
    properties 
        objective                   % objective function or struct containing function, scaling etc 
        constraints   = [];         % non-linear constraint function(s)or struct containing function(s), scaling etc 
        maps                        % mappings from/to schedule/well structures to/from scaled control-vector
        parameters                  % cell-array of parameter-objects
        bounds                      % bound constraints for controls 
        objectiveStatFun  = @mean;  % function for computing main objective based on member/realization objectives
        constraintStatFun  = @mean;  
        plotProgress  = false;  %
        solverFun               % 
        updateProblemFun  = []; % Function for updating problem structure from decision variables u
                                % of the form:  problem = updateFn(problem, u);
        fetchVariablesFun = []; % Function for fetching decision variables from problem
                                % of the form: u = fetchVariablesFun(problem) 
        
        mainDirectory  
    end
    properties (SetAccess = protected)
        %background    = false; % 	 % if true, assume run in background, disable plotting
        iterationObjectiveValues     % ResultHandler to objective for each iteration (value/gradient)
        iterationControls       % ResultHandler to "control" for each iteration
        memberObjectiveValues        % ResultHandler to objective for each ensemble member for current iteration (value/gradient)
        memberConstraintValues
        iterationConstraintValues
        realizations
        regularization 
    end
    
    methods
        function p = OptimizationProblem(samples, varargin)
            opt1 = struct('directory',   '', ...
                         'name',      'tmp');
            [opt1, rest] = merge_options(opt1, varargin{:});
            if isempty(opt1.directory)
                opt1.directory = fullfile(mrstOutputDirectory(), opt1.name);
            end
            % select first iteration directory
            firstDir = fullfile(opt1.directory, '1');
            opt = struct('objective',         [], ...
                         'constraints',       [], ...
                         'solverFun',         [], ...
                         'maps',              [], ...
                         'parameters',        [], ...
                         'bounds',            [], ...
                         'statFun',        @mean, ...
                         'updateProblemFun',  [], ...
                         'fetchVariablesFun', [], ...
                         'plotProgress',    true, ...
                         'setupType',         '', ...
                         'regularization',    []);
            [opt, rest] = merge_options(opt, rest{:});
            
            p = p@BaseEnsemble(samples, rest{:}, ...
                                 'prepareSimulation', false, ...
                                 'reset',             false, ...
                                 'directory',         firstDir, ...
                                 'name',              opt1.name, ...
                                 'solve',             []);
            p.mainDirectory = opt1.directory;
            
            p = setupProblem(p, opt);
            
            getHandler = @(nm, loc)ResultHandler('dataPrefix', nm, ...
                                                 'dataDirectory', p.mainDirectory, ...
                                                 'dataFolder', loc);
            p.memberObjectiveValues    = getHandler('objective', '1');
            p.iterationObjectiveValues = getHandler('objective', '');
            p.iterationControls        = getHandler('controls',  '');
            if ~isempty(p.constraints)
                p.memberConstraintValues    = getHandler('constraints', '1');
                p.iterationConstraintValues = getHandler('constraints', '');
            end
            
            if isempty(p.realizations)
                p.realizations = 1:p.samples.num;
            end   
        end
        
        %------------------------------------------------------------------
        %------------------------------------------------------------------
        function problem = getSampleProblem(p, seed, varargin)
            % check if problem(s) is/are explicitly given:
            if isfield(p.samples, 'getSampleProblem') || ...
                    ismethod(p.samples, 'getSampleProblem')
                problem = getSampleProblem@BaseEnsemble(p, seed, varargin{:});
            else
                % assume problem(s) are given explicitly
                problem = p.samples.problem{seed};
                problem = resetHandlers(problem, p.directory, num2str(seed));
            end
            % Setup problem for current iteration. 
            problem = p.updateProblemForCurrentIteration(problem);
        end
        
        %------------------------------------------------------------------
        function problem = updateProblemForCurrentIteration(p, problem)
            cNo = p.getCurrentControlNo();
            u   = p.unscaleVariables(p.iterationControls{cNo});
            if isempty(p.updateProblemFun)
                problem.u = u;
            else
                problem = p.updateProblemFun(problem, u);
            end
        end 
        
        %------------------------------------------------------------------
        %----------- OBJECTIVE EVALUATIONS --------------------------------
        %------------------------------------------------------------------
        
        function varargout = getScaledObjective(p, us, isConstraint)
            % Main routine used by optimizer
            if nargout < 3
                isConstraint = false;
            end
            objScaling = 1;          
            if ~isConstraint && isstruct(p.objective) && isfield(p.objective, 'scaling')
                objScaling = p.objective.scaling;
            elseif isstruct(p.constraints) && isfield(p.constraints, 'scaling')
                objScaling = p.constraints.scaling;
            end
            if nargout < 2
                varargout{1} = p.getObjective(us, isConstraint)/objScaling;
            else
                [v, dv] = p.getObjective(us, isConstraint);
                varargout{1} = v/objScaling;
                [lower, upper] = p.getControlVectorLimits();
                varargout{2} = (dv.*(upper-lower))/objScaling;
            end
            if ~isempty(p.regularization)
                [r, gr] = p.regularization(us);
                varargout{1} = varargout{1} + r;
                if nargout > 1
                    varargout{2} = varargout{2} + gr;
                end
            end
        end
        
        %------------------------------------------------------------------
        function dv = getScaledObjectiveGradient(p, control, varargin)
            % Standalone gradient routine 
            [~, dv] = p.getScaledObjective(control, varargin{:});
        end
        
        %------------------------------------------------------------------
        function varargout = getObjective(p, us, isConstraint)
            % Main unscaled routine 
            if nargin < 3
                isConstraint = false;
            end
            if ~isConstraint % default objective
                h = p.iterationObjectiveValues;
            else             % non-linear constraints
                h = p.iterationConstraintValues;
            end
            getGradient = nargout > 1;
            p    = p.processControlVector(us);
            itNo = p.getCurrentControlNo();
            isComputed = false;
            if ismember(itNo, h.getValidIds())
                obj = h{itNo};
                if ~getGradient
                    isComputed = true;
                else
                    isComputed = isfield(obj, 'gradient') && ~isempty(obj.gradient);
                    if ~isComputed
                        % objective is computed, but not gradient. Reset
                        % status:
                        failed   = p.getSimulationStatus() == -1;
                        p.simulationStatus.resetData(find(~failed));
                    end
                end
            end
            if ~isComputed
                p.computeEnsembleObjectives(getGradient, isConstraint);
                obj = p.assembleIterationObjective(isConstraint);
                h{itNo} = obj; %# ok (write file) 
            end
            if nargout >= 1
                varargout{1} = obj.value;
            end
            if nargout >= 2
                varargout{2} = obj.gradient;
            end
        end   
        
        %------------------------------------------------------------------
        function computeEnsembleObjectives(p, computeGradient, isConstraint)
            % Compute objective/constraints for each ensemble member
            if nargin < 3
                isConstraint = false;
            end
            if ~isConstraint % default objective
                [objFun, h] = deal(p.objective, p.memberObjectiveValues);
            else             % non-linear constraints
                [objFun, h] = deal(p.constraints, p.memberConstraintValues);
            end 
            if isstruct(objFun)
                objFun = objFun.function;
            end
            % set solve-function depending on what we want to compute
            p.solve = @(problem)p.solverFun(problem, objFun, h, computeGradient);
            % disregard previously failed simulations 
            failed   = p.getSimulationStatus() == -1;
            assert(~all(failed), 'All previous simulations have failed, aborting.')
            computed = false(size(failed));
            computed(h.getValidIds) = true;
            doCompute = ~computed & ~failed;
            % reset simulation status if needed
            p.simulationStatus.resetData(doCompute(p.simulationStatus.getValidIds()));
            bg = ~strcmp(p.simulationStrategy, 'serial');
            if bg
                p.prepareEnsembleSimulation();
            end
            p.simulateEnsembleMembers('range', reshape(p.realizations(doCompute), 1, []));
            if bg
                % wait for computations to complete
                isDone = false;
                while ~isDone
                    isDone = nnz(~failed) == numel(p.memberObjectiveValues.getValidIds);
                    pause(1);
                    % some plotting?
                end
            end
        end
        
        % -----------------------------------------------------------------
        function obj = assembleIterationObjective(p, isConstraint)
            % compute mean (or other statistic) objective for current
            % iteration
            if nargin < 2
                isConstraint = false;
            end
            if ~isConstraint
                h = p.memberObjectiveValues;
                statfun = p.objectiveStatFun;
            else % non-linear constraint
                h = p.memberConstraintValues;
                statfun = p.constraintStatFun;
            end
            ix = find(p.getSimulationStatus == 1);
            % do simple check
            tmp = false(p.samples.num, 1);
            tmp(h.getValidIds()) = true;
            assert(all(tmp(ix)==1), 'Unexpected non-valid handler ids')
            %
            objList = h(ix);
            vals = cellfun(@(x)x.value, objList);
            okVals = isfinite(vals);
            if any(okVals)
                obj.value = statfun(vals(okVals));
            else
                warning('Not able to produce objective value for iteration %d', p.getCurrentControlNo);
                obj.value = nan;
            end
            if ~isempty(objList) && isfield(objList{1}, 'gradient')
                grads = applyFunction(@(x)x.gradient, objList);
                okGrads = cellfun(@(x)~isempty(x) && all(isfinite(x)), grads);
                if any(okGrads)
                    obj.gradient = statfun(horzcat(grads{:}),2);
                end
            end         
        end
        
        %------------------------------------------------------------------
        %------- UTILITIES  -----------------------------------------------
        %------------------------------------------------------------------

        function cNo = getCurrentControlNo(p)
            [~, s] = fileparts(p.directory);
            cNo    = str2double(s);
            assert(~isnan(cNo), 'Unexpected directory name: %s\n', p.directory)
        end
        
        %------------------------------------------------------------------
        function us = scaleVariables(p, u)
            [lower, upper] = p.getControlVectorLimits();
            us = (u-lower)./(upper-lower);
        end
        
        %------------------------------------------------------------------
        function u = unscaleVariables(p, us)
            [lower, upper] = p.getControlVectorLimits();
            u = lower + us.*(upper-lower);
        end
        
        %------------------------------------------------------------------
        function [lower, upper] = getControlVectorLimits(p)
            if ~isempty(p.bounds)
                [lower, upper] = deal(p.bounds(:,1), p.bounds(:,2));
            else
                [lower, upper] = deal(0, 1); % no effect
            end
        end
          
        %------------------------------------------------------------------
        function p = processControlVector(p, us)
            % Identify if a control-vector already exist and if so switch
            % to that folder. If new, then setup new folder (and switch) 
            controlNo = [];
            vc = p.iterationControls.getValidIds();
            if any(vc)
                ix = flipud(vc(:));
                for k = 1:numel(ix)
                    if norm(single(us)-single(p.iterationControls{ix(k)})) == 0
                        controlNo = ix(k);
                        break;
                    end
                end
            end
            if isempty(controlNo)
                if isempty(vc)
                    controlNo = 1;
                else
                    controlNo = max(vc) +1;
                end
                p.iterationControls{controlNo} = {us};
            end
            folderNm = fullfile(p.mainDirectory, num2str(controlNo));
            if ~exist(folderNm, 'dir')
                mkdir(folderNm);
            end
            % switch to directory
            p = switchDirectory(p, controlNo);
        end
        
        %------------------------------------------------------------------
        function p = switchDirectory(p, cNo)
            assert(isnumeric(cNo), 'Directory name is expected to be of numeric type')
            dd = num2str(cNo);
            folderNm = fullfile(p.mainDirectory, dd);
            if exist(folderNm, 'dir')
                p.directory = folderNm;
                p.simulationStatus.dataFolder = dd;
                p.memberObjectiveValues.dataFolder = dd;
                if ~isempty(p.memberConstraintValues)
                    p.memberConstraintValues.dataFolder = dd;
                end
            else
                error('Non-existent directory: %s\n', folderNm);
            end
        end
        
        % -----------------------------------------------------------------
        function flag = reset(p, varargin)
            opt = struct('prompt', true, 'prepareSimulation', true);
            opt = merge_options(opt, varargin{:});
            if ~isempty(p.iterationControls.getValidIds())
                if opt.prompt
                    prompt = sprintf(['Delete all data for %s? (sample '    , ...
                        'data will not be deleted) y/n [n]: '], p.name);
                    if ~strcmpi(input(prompt, 's'), 'y')
                        fprintf('Ok, will not remove files.\n');
                        flag = false;
                        return
                    end
                end
                ids = p.iterationControls.getValidIds();
                for cNo = 1:max(ids)
                    p = p.switchDirectory(cNo);
                    p.memberObjectiveValues.resetData();
                    flag = reset@BaseEnsemble(p, 'prompt', false, 'prepareSimulation', false);
                end
                p.directory = fullfile(p.mainDirectory, '1');
                p.iterationObjectiveValues.resetData();
                p.iterationControls.resetData();
                if opt.prepareSimulation
                    p.prepareEnsembleSimulation();
                end
            end
        end 
        
        %------------------------------------------------------------------
        %------- OPTIMIZATION  --------------------------------------------
        %------------------------------------------------------------------
        function [u, extra] = maximizeObjective(p, varargin)
            [u, extra] = optimize(p, varargin{:}, 'maximize', true);
        end
        
        function [u, extra] = minimizeObjective(p, varargin)
            [u, extra] = optimize(p, varargin{:}, 'maximize', false);
        end
        
        function [u, extra] = optimize(p, varargin)
            if mod(nargin, 2) == 0 % first argument is initial guess
                initialGuess = varargin{1};
                varargin     = varargin(2:end);
                % initialGuess is either assumed to be an unscaled vector
                % of decision varables, or a problem structure
                if ~isa(initialGuess, 'double')
                    try
                        initialGuess = p.fetchVariablesFun(initialGuess);
                    catch
                        error('Not able to fetch decision variables from initialGuess');
                    end
                end
                us0 = p.scaleVariables(initialGuess);
            else
                % assume valid initial exists on file
                ids = p.iterationControls.getValidIds();
                assert(~isempty(ids) && ids(1)==1, ...
                    'Not able to fetch initial guess');
                us0 = p.iterationControls{1};
            end
            assert(all(-sqrt(eps)<us0 & us0 < 1+sqrt(eps)), ...
                'Initial guess is not within bounds');
            opt = struct('optimizer',      'default', ...
                         'background',     false, ...
                         'maximize',       true,  ...
                         'regularizationParam', 0);
            [opt, optimopts] = merge_options(opt, varargin{:});
            if opt.regularizationParam > 0
               p.regularization = getRegularizationFun(us0, opt.regularizationParam);
            else
               p.regularization = [];
            end
            if ~opt.background
                switch lower(opt.optimizer)
                    case 'default'
                        [~, us, hist] = unitBoxBFGS(us0, @p.getScaledObjective, optimopts{:}, 'maximize', opt.maximize);
                        if nargout == 2
                            extra = hist;
                        end
                        
                    case 'ipopt'
                        nu = numel(us0);
                        objSign = 1;
                        if opt.maximize
                            objSign = -1;
                        end
                        funcs = struct('objective',    @(u)objSign*p.getScaledObjective(u'), ...
                                       'gradient',     @(u)objSign*p.getScaledObjectiveGradient(u').');
                        if ~isempty(p.constraints)
                            funcs.constraints = @(u)p.getScaledObjective(u, p.constraint, true);
                            funcs.jacobian    = @(u)p.getScaledObjectiveGradient(u, p.constraint, true).';
                        end
                        options = struct('lb', zeros(1, nu), 'ub', ones(1, nu));
                        
                        ipopt_opts = struct('hessian_approximation', 'limited-memory', ...
                                            'tol',                   1e-4, ...
                                            'max_iter',              25);
                        ipopt_opts = merge_options(ipopt_opts, optimopts{:});
                        options.ipopt = ipopt_opts;
                        [us, info] = ipopt(us0.',funcs , options);
                        us = us(:);
                        if nargout == 2
                            extra = info;
                        end
                    otherwise
                        error('Unsupported optimizer: %s', opt.optimizer)
                end
                u = p.unscaleVariables(us);
            else
                p.optimizeBackground(us0, [{'optimizer', opt.optimizer}, optimopts]);
                [u, extra] = deal(true, 'Started optimization in seperate session');
            end
            
        end
        
        %------------------------------------------------------------------
        function [] = optimizeBackground(p, u0, varargin)
            if ~ismember(1, p.iterationControls.getValidIds())
                p.iterationControls{1} = {u0};
            end
            optimOpts = [varargin{:}, {'plotEvolution', false}];
            % We don't know what optimOpts are, save in case wrapper can't
            % handle them:
            optname = fullfile(p.mainDirectory, 'optimOpts.mat');
            save(optname, 'optimOpts');
            % save problem
            fname = fullfile(p.mainDirectory, 'optimizationProblem.mat');
            save(fname, 'p');
            % Run optimization in seperate matlab session
            fprintf('Starting optimization of problem %s in background\n', p.name);
            % call with correct class-name
            evalFunWrapper(@optimizeStandalone, {fname, optname}, 'exitWhenDone', false);
        end
        
        
        %------------------------------------------------------------------
        %------- plotting functionalites ----------------------------------
        %------------------------------------------------------------------
        function h = plotObjectiveValues(p, includeMembers)
            if nargin < 2
                includeMembers = true;
            end
            controlIds = p.iterationControls.getValidIds();
            nc = max(controlIds);
            if nc > 0
                vals = nan(1, nc);
                ids  = p.iterationObjectiveValues.getValidIds;
                vals(ids) = cellfun(@(x)x.value, p.iterationObjectiveValues(ids));
                h(1) = plot(1:nc, vals, '-ob', 'LineWidth', 3);
                if includeMembers
                    hold on
                    curIt = p.getCurrentControlNo;
                    V = nan(nc, p.samples.num);
                    for k = 1:nc
                        p.memberObjectiveValues.dataFolder = num2str(k);
                        ids = p.memberObjectiveValues.getValidIds;
                        V(k, ids) = cellfun(@(x)x.value, p.memberObjectiveValues(ids));
                    end
                    p.memberObjectiveValues.dataFolder = num2str(curIt);
                    h(1+(1:p.samples.num)) = plot((1:nc)', V, '-o');
                    nms = applyFunction(@(n)sprintf('realization %d', n), ...
                        p.realizations);
                    legend(h, 'main', nms{:});
                    title('Objective values')
                    xlabel('Evaluation number')
                end
            end
        end
        
        %------------------------------------------------------------------
        function [] = plotWellSols(p, problem, controlNo, memberNo)
            if nargin < 4
                memberNo = p.realizations;
            end
            [nc,nm] = deal(numel(controlNo),numel(memberNo));
            [wss, tms, nms] = deal(cell(1, nc*nm));
            cnt     = 0;
            tm = problem.SimulatorSetup.schedule.step.val;
            %cc = p.getCurrentControlNo();
            for kc = nc
                for km = 1:nm
                    cnt = cnt +1;
                    [cno, memno] = deal(controlNo(kc), memberNo(km));
                    %p.switchDirectory(cno);
                    wsh = getHandler(p.mainDirectory, 'wellSols', cno, memno, false);
                    if ~isempty(wsh)
                        ids = wsh.getValidIds;
                        wss{cnt} = wsh(ids);
                        wss{cnt} = wss{cnt}(:);
                        tms{cnt} = tm(ids);
                        nms{cnt} = sprintf('Schedule %d, realization %d', cno, memno);
                    end
                end
            end
            ix = ~cellfun(@isempty, wss);
            if ~any(ix)
                warning('No well-solutions found for the requested controls/realizations\n')
            else
                if ~all(ix)
                    fprintf('Found well-solutions from only some of the requested controls/realizations\n')
                end
                plotWellSols(wss(ix), tms(ix), 'Datasetnames', nms(ix));
            end
        end
        
        function problem = resetHandlers(p, problem, directory,folder)
            problem = resetHandlers(problem, directory, folder);
        end
    end
end

% -----------------------------------------------------------------
% Helper functions for well-controls (non-empty 'maps'-prop)
% -----------------------------------------------------------------

function problem = updateWellControls(problem, u, maps)
schedule = problem.SimulatorSetup.schedule;
nc = numel(maps.type);
for k = 1:nc
    [sno, wno, tp] = deal(maps.stepNo(k), maps.wellNo(k), maps.type{k});
    if maps.isTarget(k)
        assert(strcmp(tp, schedule.control(sno).W(wno).type), ...
            'Schedule contains non-expected well type');
        schedule.control(sno).W(wno).val = u(k);
        % also update if part of lims
        if isfield(schedule.control(sno).W(wno).lims, tp)
            schedule.control(sno).W(wno).lims.(tp) = u(k);
        end
    else
        schedule.control(sno).W(wno).lims.(tp) = u(k);
    end
end
problem.SimulatorSetup.schedule = schedule;
if isfield(schedule.control(1).W(1), 'posControl') && nc < numel(u)
    % trajectory controls
    W = schedule.control(1).W;
    pcix = arrayfun(@(w)~isempty(w.posControl), W);
    pc = {W(pcix).posControl};
    ix = nc;
    for k = 1:numel(pc)
        nParam = pc{k}.parameters.nParam;
        pc{k} =  pc{k}.control2param(u(ix + (1:nParam)));
        ix = ix + nParam;
    end
end
end
% -----------------------------------------------------------------

function u = fetchWellControlValues(problem, maps)
u = zeros(numel(maps.type), 1);
schedule = problem.SimulatorSetup.schedule;
for k = 1:numel(u)
    [sno, wno, tp] = deal(maps.stepNo(k), maps.wellNo(k), maps.type{k});
    if maps.isTarget(k)
        assert(strcmp(tp, schedule.control(sno).W(wno).type), ...
            'Schedule contains non-expected well type');
        u(k) = schedule.control(sno).W(wno).val;
    else
        u(k) = schedule.control(sno).W(wno).lims.(tp);
    end
end

if isfield(schedule.control(1).W(1), 'posControl')
    % position is same over all steps
    W = schedule.control(1).W;
    ix = arrayfun(@(w)~isempty(w.posControl), W);
    pc = {W(ix).posControl};
    up = cellfun(@(x)x.param2control(), pc, 'UniformOutput', false);
    u = [u; vertcat(up{:})];
end
end

% -------------------------------------------------------------------------
% Helper functions for parameter-problems (non-empty 'parameters'-prop)
% -------------------------------------------------------------------------

function problem = updateModelParameters(problem, u, params)
nparam = cellfun(@(x)x.nParam, params);
u = mat2cell(u, nparam, 1);
% Create new setup, and set parameter values
setup = problem.SimulatorSetup;
% make sure to recreate discretizations
setup.model.FlowDiscretization = [];
setup.model.FlowPropertyFunctions = [];
for k = 1:numel(params)
    tmp   = params{k}.unscale(u{k});
    setup = params{k}.setParameter(setup, tmp);
end
problem.SimulatorSetup = setup;
end
% -------------------------------------------------------------------------

function u = fetchParameterValues(problem, params)
u = getScaledParameterVector(problem.SimulatorSetup, params);
end


% -------------------------------------------------------------------------
% Setup of default simulation/diagnsotics problems
% -------------------------------------------------------------------------
function p = setupProblem(p, opt)
props = intersect(properties(p), fieldnames(opt));
for k = 1:numel(props)
    p.(props{k}) = opt.(props{k});
end

if ~isempty(opt.setupType)
    % setup simulation/diagnostics problem
    isWellType  = ~isempty(p.maps);
    isParamType = ~isempty(p.parameters);
    assert(isWellType || isParamType, ...
        'Default setup requires either input of ''maps'' or ''parameters''');
    if isWellType && isParamType
        warning('Can''t work with both well controls and parameters. %s', ...
                'Parameters will be disregarded');
        p.parameters = [];
    end
    if isempty(p.solverFun)
        if strcmp(opt.setupType, 'simulation')
            p.solverFun = @(problem, obj, h, computeGradient) ...
                simulationSolverFun(problem, obj, 'objectiveHandler', h, ...
                'computeGradient', computeGradient, ...
                'parameters', p.parameters, 'maps', p.maps);
        elseif strcmp(opt.setupType, 'diagnostics')
            p.solverFun = @(problem, obj, h, computeGradient) ...
                diagnosticsSolverFun(problem, obj, 'objectiveHandler', h, ...
                'computeGradient', computeGradient, ...
                'parameters', p.parameters);
        else
            error('Unknown setup-type: %s', opt.setupType);
        end
    end
    if isempty(p.updateProblemFun)
        if isWellType    
            p.updateProblemFun = ...
                @(problem,u)updateWellControls(problem, u, p.maps);
        else % parameter type
            p.updateProblemFun = ...
                @(problem,u)updateModelParameters(problem, u, p.parameters);
        end
    end
    if isempty(p.fetchVariablesFun)
        if isWellType 
            p.fetchVariablesFun = ...
                @(problem)fetchWellControlValues(problem, p.maps);
        else% parameter problem
            p.fetchVariablesFun = ...
                @(problem)fetchParameterValues(problem, p.parameters);
        end
    end
    % deal with bounds
    if isWellType
        if ~isempty(p.bounds)
            warning('Disregarding input: ''bounds'' (using maps.bounds)');
        end
        p.bounds = p.maps.bounds;
    end
    if isParamType
        if ~isempty(p.bounds)
            warning('Disregarding input: ''bounds'' (using parameter.boxLims)');
        end
        p.bounds = [];
    end
end


end

% -------------------------------------------------------------------------
% Handler utilities
% -------------------------------------------------------------------------

function h = getHandler(maindd, prefix, controlNo, memberNo, showWarning)
if nargin < 5
    showWarning = true;
end
pth = fullfile(maindd, num2str(controlNo), num2str(memberNo));
if isfolder(pth)
    [dd, fld] = fileparts(pth);
    h = ResultHandler('dataPrefix', prefix, ...
        'writeToDisk', true,...
        'dataDirectory', dd, ...
        'dataFolder', fld, ...
        'cleardir', false);
else
    h = [];
    if showWarning
        warning('Non-existing directory: %s', pth)
    end
end
end

function problem = resetHandlers(problem, directory, folder)
if isfield(problem, 'OutputHandlers')
    fn = fieldnames(problem.OutputHandlers);
    for k = 1:numel(fn)
        h = problem.OutputHandlers.(fn{k});
        if isa(h, 'ResultHandler')
            if ~exist(fullfile(directory, folder), 'dir')
                mkdir(fullfile(directory, folder));
            end
            h.dataDirectory = directory;
            h.dataFolder    = folder;
            problem.OutputHandlers.(fn{k}) = h;
        end
    end
end
end
    
% -------------------------------------------------------------------------
% Regularization
% -------------------------------------------------------------------------

function fn = getRegularizationFun(us0, alpha)
fn = @(us)regularizationFn(us, us0, alpha);
end

function varargout = regularizationFn(us, us0, alpha)
varargout{1} = alpha*sum( (us-us0).^2 );
if nargout > 1
    varargout{2} = 2*alpha*(us-us0);
end
end