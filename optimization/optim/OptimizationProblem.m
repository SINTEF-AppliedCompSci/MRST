classdef OptimizationProblem < handle
    properties 
        name                    % problem name
        objective               % function handle to objective function
        objectiveOpts  = {};    % optional input to objective
        realizations  =  1;     % realization numbers to consider
        initialGuess            % initial guess for controls (for simulation-problems: schedule) 
        setupFn                 % function handle to setup-function (for simulation problems: setup(n) -> problem)
        evalFn                  % function handle to objective/gradient-function (default assumption: simulations)
        evalType =   '';        % type of evaluation (default 'simulation')          
        maps                    % mappings from/to schedule/well structures to/from scaled control-vector
        bounds                  % bound constraints for controls
        constraint              % optional function handle to a (lumped) constraint function (for use with external optimizers) 
        outputPath              % optimization problem main output-path
        nWorkers = 0;           % number of external matlab-processes used for computations (zero implies only current)
        optimizer               % name (string) of optimizer. Default is the mrst-option unitBoxBFGS
        objectiveScaling = 1;   % objective scaling. A good guess of objective magnitude helps the optimizer choose a reasonnable initial step-length.
        plotProgress  = false;   %
        background    = false; % if true, assume run in background, disable plotting
    end
    properties (SetAccess = private)
        evaluationFolders  = {};     % list of all created evaluation-output-folders     
        controlList        = [];     % list of all scaled controls
    end
    
    methods
        function p = OptimizationProblem(name, varargin)
            p.name = name;
            % check for non-default output-path
            ix = find(strcmp('outputPath', varargin));
            if ~isempty(ix)
                p.outputPath = varargin{ix+1};
            else
                p.outputPath = fullfile(mrstOutputDirectory(), name);
            end
            
            % check if problem allready has been created 
            probFile = fullfile(p.outputPath, 'problemSetup.mat');
            if isfile(probFile)
                fprintf('Loading previously created optimization problem setup')
                tmp = load(probFile);
                flds = fieldnames(tmp.p);
                for k = 1:numel(flds)
                    p.(flds{k}) = tmp.p.(flds{k});
                end
                clear tmp
            end
           
            [p, optimOpts] = merge_options(p, varargin{:});
            % check what type of evaluations we are dealing with
            if isempty(p.evalType)
                if isa(p.objective, 'DiagnosticsObjective')
                    p.evalType = 'diagnotics';
                else
                    p.evalType = 'simulation';
                    p.maps   = setupSimulationControlMappings(p.initialGuess, p.bounds); 
                    p.evalFn = @computeEnsembleSimulationObjective;
                end
            end
            
            p = p.checkProblemSetup();
            
            % setup folder ------------------------------------------------      
            if ~isfolder(p.outputPath)
                mkdir(p.outputPath);
            end
            
            % save setup
            if ~p.background
                save(fullfile(p.outputPath, 'problemSetup.mat'), 'p');
            end
            
            % get previously computed output
            p = p.refreshControlList();
            
            % if remaining options are non-empty, optimization is expected to be
            % run straight away (possibly called from different session)
            if ~isempty(optimOpts)
                opt = struct('runOptimization', false);
                [opt, remaining] = merge_options(opt, optimOpts{:});
                if opt.runOptimization
                    p.optimize(remaining{:});
                else
                    error('Unexpected input');
                end
            end
        end
        
        %------------------------------------------------------------------
        
        function str = evalOutputPath(p, simNo)
            str = fullfile(p.outputPath, num2str(simNo));
        end
        
        %------------------------------------------------------------------
        
        function str = memberOutputPath(p, simNo, memNo)
            str = fullfile(p.evalOutputPath(simNo), num2str(memNo));
        end
        
        %------------------------------------------------------------------
        
        function value = getMemberObjective(p, obj, memberNo, varargin)
            % get/compute objective for given control for realization memberNo
            [~, ~, curpath] = p.handleControlInfo(memberNo, varargin{:});
            fileNm = fullfile(curpath, [func2str(p.objective), '.mat']);
            if isfile(fileNm) 
                tmp = load(fileNm);
                value = tmp.value;
            else
                objOpt = p.getObjectiveOptions(obj);
                value = p.evalFn(p.setupFn, 'outputPath', fileparts(curpath), ...
                                  'objective', obj, 'computeGradient', false, ...
                                  'memberIx', memberNo, objOpt{:});
                save(fileNm, 'value');
            end
        end
        
        %------------------------------------------------------------------
        
        function value = getObjective(p, obj, varargin)
            % get/compute (expected) objective for given control
            [control, controlNo, curpath] = p.handleControlInfo([], varargin{:});
            fileNm = fullfile(curpath, [func2str(p.objective), '.mat']);
            if isfile(fileNm)
                tmp = load(fileNm);
                value = tmp.value;
            else
                if p.nWorkers > 0
                    % handle computations in seperate matlab sessions
                    objOpt = p.getObjectiveOptions(obj);
                    args = [{'outputPath',      curpath, ...
                            'objective',       obj, ...
                            'computeGradient', true}, objOpt];
                    p.launchBackgroundSessions(p.evalFn, args, curpath)
                end
                % compute or read individual values
                [value, numOK] = deal(0);
                for k = 1:numel(p.realizations)
                    memberNo = p.realizations(k);
                    vi = p.getMemberObjective(obj, memberNo, 'control', control, 'controlNo', controlNo);
                    if ~isempty(vi)
                        numOK = numOK +1;
                        value = value + vi;
                    end
                end
                value = value/numOK;
                save(fileNm, 'value');
                if p.plotProgress
                    ax = gca; cla(ax);
                    p.plotObjectiveValues();
                end
            end
        end
        
        %------------------------------------------------------------------
        
        function gradient = getMemberGradient(p, obj, memberNo, varargin)
            % get/compute gradient of given control for realization memberNo
            [~, ~, curpath] = p.handleControlInfo(memberNo, varargin{:});
            fileNm = fullfile(curpath, ['gradient_', func2str(p.objective), '.mat']);
            if isfile(fileNm)
                tmp = load(fileNm);
                gradient = tmp.gradient;
            else
                if ~isfolder(curpath)
                    mkdir(curpath);
                end
                objOpt = p.getObjectiveOptions(obj);
                [~, gradient] = p.evalFn(p.setupFn, 'outputPath', fileparts(curpath), ...
                                  'objective', obj, 'computeGradient', true, ...
                                  'memberIx', memberNo, objOpt{:});
                save(fileNm, 'gradient');
            end
        end
        
        %------------------------------------------------------------------
        
        function gradient = getGradient(p, obj, varargin)
            % get/compute gradient of (expected) objective for given control
            [control, controlNo, curpath] = p.handleControlInfo([], varargin{:});
            fileNm = fullfile(curpath, ['gradient_', func2str(p.objective), '.mat']);
            if isfile(fileNm)
                tmp = load(fileNm);
                gradient = tmp.gradient;
            else
                if p.nWorkers > 0
                    % handle computations in seperate matlab sessions
                    objOpt = p.getObjectiveOptions(obj);
                    args = [{'outputPath',      curpath, ...
                            'objective',       obj, ...
                            'computeGradient', true}, objOpt];
                    p.launchBackgroundSessions(p.evalFn, args, curpath)
                end
                % compute or read individual values
                [gradient, numOK] = deal([], 0);
                for k = 1:numel(p.realizations)
                    memberNo = p.realizations(k);
                    gi = p.getMemberGradient(obj, memberNo, 'control', control, 'controlNo', controlNo);
                    if ~isempty(gi)
                        numOK = numOK + 1;
                        if isempty(gradient)
                            gradient = gi;
                        else
                            if isa(gradient, 'double')
                                gradient = gradient + gi;
                            else
                                flds = fieldnames(gradient);
                                for nk = 1:numel(flds)
                                    gradient.(flds{nk}) = gradient.(flds{nk}) + gi.(flds{nk});
                                end
                            end
                        end
                    end
                end
                if isa(gradient, 'double')
                    gradient = gradient/numOK;
                else
                    flds = fieldnames(gradient);
                    for nk = 1:numel(flds)
                        gradient.(flds{nk}) = gradient.(flds{nk})/numOK;
                    end
                end
                save(fileNm, 'gradient')
            end
        end
        
        %------------------------------------------------------------------
        
        function varargout = getScaledObjective(p, control, obj)
            % get/compute objective (and possibly gradient) of scaled
            % problem. This is the main evalution function that should be 
            % handed to the optimizer
            if nargin < 3
                obj   = p.objective;
            end
            scale = 1;
            if isequal(obj, p.objective)
                scale = p.objectiveScaling;
            end
            varargout{1} = p.getObjective(obj, 'control', control)/scale;
           
            if nargout > 1
                varargout{2} = p.getScaledGradient(control, obj);
            end
        end
        
        %------------------------------------------------------------------        
        
        function g = getScaledMemberGradient(p, obj, memberNo, varargin)
            % get/compute gradient of scaled problem for realization memberNo 
            [control, controlNo] = p.handleControlInfo(memberNo, varargin{:});
            gradient = p.getMemberGradient(obj, memberNo, 'control', control, 'controlNo', controlNo);
            if strcmp(p.evalType, 'simulation')
                g = mapGradient(gradient, p.maps);
            end
            scale = 1;
            if isequal(obj, p.objective)
                scale = p.objectiveScaling;
            end
            g = g/scale;
        end
        
        %------------------------------------------------------------------
        
        function g = getScaledGradient(p, control, obj)
            % get/compute gradient of the scaled problem. This function should 
            % be handed to optimizers requiring a seperate gradient-function  
            if nargin < 3
                obj   = p.objective;
            end
            [control, controlNo] = p.handleControlInfo([], 'control', control);
            % We the call gradient function just to make sure everything 
            % gets computed (the rigt place)
            getGradient(p, obj, 'control', control, 'controlNo', controlNo);
            % We need to process and scale all member gradients
            % compute or read individual values
            [g, numOK] = deal([], 0);
            for k = 1:numel(p.realizations)
                memberNo = p.realizations(k);
                gi = p.getScaledMemberGradient(obj, memberNo, 'control', control, 'controlNo', controlNo);
                if ~isempty(gi)
                    numOK = numOK + 1;
                    if isempty(g)
                        g = gi;
                    else
                        if isa(g, 'double')
                            g = g + gi;
                        else
                            flds = fieldnames(g);
                            for nk = 1:numel(flds)
                                g.(flds{nk}) = g.(flds{nk}) + gi.(flds{nk});
                            end
                        end
                    end
                end
            end
            if isa(g, 'double')
                g = g/numOK;
            else
                flds = fieldnames(g);
                for nk = 1:numel(flds)
                    g.(flds{nk}) = g.(flds{nk})/numOK;
                end
            end
        end
        
        %------------------------------------------------------------------
        
        function [schedule, controlNo, extra] = optimize(p, varargin)
            % run optimization
            opt = struct('optimizer', '');
            [opt, optimopts] = merge_options(opt, varargin{:});
            if isempty(opt.optimizer)
                opt.optimizer = 'unitboxBFGS';
            end
            if strcmp(p.evalType, 'simulation')
                u0 = p.schedule2control(p.initialGuess);
            else
                u0 = p.initialGuess;
            end
            switch lower(opt.optimizer)
                case 'unitboxbfgs'
                    [~, u, hist] = unitBoxBFGS(u0, @p.getScaledObjective, optimopts{:}, 'plotEvolution', ~p.background);
                    if nargout == 3
                        extra = hist;
                    end
                    
                case 'ipopt'
                    nu = numel(u0);
                    func = struct('objective',    @(u)-p.getScaledObjective(u), ...
                                  'gradient',     @(u)-p.getScaledGradient(u).');
                    if ~isempty(p.constraint)
                        func.constraints = @(u)p.getScaledObjective(u, p.constraint);
                        func.jacobian    = @(u)p.getScaledGradient(u, p.constraint).';
                    end
                    options = struct('lb', zeros(1, nu), 'lu', ones(1, nu));

                    ipopt_opts = struct('hessian_approximation', 'limited-memory', ...
                                        'tol',                   1e-4, ...
                                        'max_iter',              25);
                    ipopt_opts = merge_options(ipopt_opts, optimopts{:});
                    options.ipopt = ipopt_opts;
                    [u, info] = ipopt(u0.',func , options);
                    u = u(:);
                    if nargout == 3
                        extra = info;
                    end
                otherwise
                    error('Unsupported optimizer: %s', opt.optimizer)
            end
            schedule = p.control2schedule(u);
            [~, controlNo] = p.handleControlInfo([], 'control', u);
        end
        
        %------------------------------------------------------------------
        function [] = optimizeBackground(p, varargin)
            % Run optimization in seperate matlab session
            args = [{p.name}, {'runOptimization', true, 'background', true}, varargin];
            fprintf('Starting optimization of problem %s in background\n', p.name);
            %evalFunWrapper(@OptimizationProblem, args, 'exitWhenDone', false, 'pathList', {pwd}); 
            evalFunWrapper(@OptimizationProblem, args); 
        end
        
       
        %------------------------------------------------------------------
        
        function launchBackgroundSessions(p, fn, args, pth)
            % lauch evaluations in a number of background sessions accoring to nWorkers 
            assert(p.nWorkers > 0); 
            nComput = numel(p.realizations);
            binSize = ceil(nComput/p.nWorkers);
            progressFiles = {}; 
            for k = 1:p.nWorkers
                bin = (1+(k-1)*binSize) : min(k*binSize,nComput);
                if ~isempty(bin)
                    curargs = [{p.setupFn}, args, {'memberIx', p.realizations(bin)}];
                    progFileNm = fullfile(pth, ['progress', num2str(k), '.mat']);
                    progressFiles = [progressFiles, {progFileNm}];  %#ok
                    fprintf('Starting new Matlab session ...\n')
                    evalFunWrapper(fn, curargs, 'exitWhenDone', true, 'progressFileNm', progFileNm);
                end
            end
            % wait until done
            fprintf('Waiting for background processes to finnish ...')
            while any(cellfun(@isfile, progressFiles))
                pause(1)
            end
            
        end
        
        %------------------------------------------------------------------
        
        function h = plotObjectiveValues(p, obj)
            if nargin < 2
                obj = p.objective;
            end
            fileNm = [func2str(obj), '.mat'];
            [mv, v] = p.loadAllFiles(fileNm, 'value');
            [mv, v] = deal(cell2mat(mv), cell2mat(v));
            nr = numel(p.realizations);
            if ~isempty(mv)
                hold on
                h(1) = plot(1:numel(mv), mv', '-o', 'LineWidth', 2);
                if nr > 1
                    h(2:(nr+1)) = plot(1:numel(mv), v', '--o');
                    nms = arrayfun(@(n)sprintf('realization %d', n), p.realizations, ...
                        'UniformOutput', false);
                    legend(h, 'mean', nms{:});
                end
                title('Objective values')
                ylabel(func2str(obj));
                xlabel('Evaluation number')
            else
                for k =1:nr+1
                    hold on
                    h(k) = plot(nan); %#ok
                end
            end
        end
        
       %-------------------------------------------------------------------
       
       function [] = monitorProgress(p, refreshRate)
           if nargin < 2
               refreshRate = 2;
           end
           h = p.plotObjectiveValues();
           
           fprintf('Updating plot every %d seconds. Close figure to quit plotting ', round(refreshRate));
           cont = true;
           cnt = 0;
           while cont
               cnt = cnt +1;
               p = p.refreshControlList();
               fileNm = [func2str(p.objective), '.mat'];
               [mv, v] = p.loadAllFiles(fileNm, 'value', false);
               [mv, v] = deal(cell2mat(mv), cell2mat(v));
               if ~isempty(mv)
                   h(1).XData = (1:numel(mv));
                   h(1).YData = mv;
                   if numel(p.realizations) > 1
                       for k = 1:size(v,1)
                           h(k+1).XData = (1:numel(v(k,:)));
                           h(k+1).YData = v(k,:);
                       end
                   end
               end
               fprintf(repmat('*', [1, mod(cnt,4)]));
               pause(refreshRate)
               fprintf(repmat('\b', [1, mod(cnt,4)]));
               cont = isvalid(h(1));
           end
           fprintf('\n')
       end
                     
       %-------------------------------------------------------------------
        
        function [] = plotWellSols(p, controlNo, memberNo)
            if nargin < 3
                memberNo = p.realizations;
            end
            [nc,nm] = deal(numel(controlNo),numel(memberNo));
            [wss, tms, nms] = deal(cell(1, nc*nm));
            cnt     = 0;
            for kc = nc
                % load schedule
                shfl = fullfile(p.evalOutputPath(controlNo), 'schedule.mat');
                if isfile(shfl)
                    tmp  = load(shfl);
                    tm   = tmp.schedule.step.val;
                    for km = 1:nm
                        cnt = cnt +1;
                        [cno, memno] = deal(controlNo(kc), memberNo(km));
                        wsh = p.getHandler('wellSols', cno, memno, false);
                        if ~isempty(wsh)
                            ids = wsh.getValidIds;
                            wss{cnt} = wsh(ids);
                            wss{cnt} = wss{cnt}(:); 
                            tms{cnt} = tm(ids);
                            nms{cnt} = sprintf('Schedule %d, realization %d', cno, memno);
                        end
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
            
        %------------------------------------------------------------------
        
        function h = getHandler(p, prefix, controlNo, memberNo, showWarning)
            if nargin < 5
                showWarning = true;
            end
            pth = p.memberOutputPath(controlNo, memberNo);
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
        
        %------------------------------------------------------------------
        
        function [control, controlNo, curpath] = handleControlInfo(p, memberNo, varargin)
            % Handle control-info. Based on input, give control-vector, control-number
            % and corresponding output-path. For new control-vector,
            % output-path will be created.
            opt = struct('control', [], 'controlNo', []);
            opt = merge_options(opt, varargin{:});
            [control, controlNo] = deal(opt.control, opt.controlNo);
            if isempty(controlNo) % we need to check exisiting controls
                assert(~isempty(control), 'Need either control or controlNo ...')
                % search for mathcing previous controls (use single precision)
                for k = numel(p.controlList):-1:1
                    if norm( single(control) - single(p.controlList{k})) == 0
                        controlNo = k;
                        break;
                    end
                end
            end
            
            % if no controlNo is found, we have a new control, record accordingly
            if isempty(controlNo) 
                controlNo   = numel(p.controlList) +1;
                controlFolder = p.evalOutputPath(controlNo);
                if ~isfolder(controlFolder)
                    mkdir(controlFolder);
                end
                % save control and add to control-list
                save(fullfile(controlFolder, 'control.mat'), 'control');
                p.controlList = [p.controlList, {control}];
                % also save make and save well/schedule
                if strcmp(p.evalType, 'simulation')
                    schedule = p.control2schedule(control);
                    save(fullfile(controlFolder, 'schedule.mat'), 'schedule');
                elseif strcmp(p.evalType, 'diagnostics')
                    W = p.control2well(control);
                    save(fullfile(controlFolder, 'W.mat'), 'W');
                end
            end
            % if schedule is empty, attempt to load it
            controlFolder = p.evalOutputPath(controlNo);
            if isempty(control)
                controlFileNm = fullfile(controlFolder, 'control.mat');
                assert(isfile(controlFileNm), 'Found no control file for the supplied scheduleNo');
                opt = load(controlFileNm);
                control = opt.control;
            end
            
            % finally find relevant folder-name
            if isempty(memberNo)
                curpath = controlFolder;
            else
                curpath = p.memberOutputPath(controlNo, memberNo);
            end
        end
        
        %------------------------------------------------------------------
        function p = checkProblemSetup(p)
            p.evalFn    = checkFunction(p.evalFn);
            p.objective = checkFunction(p.objective);
            p.setupFn   = checkFunction(p.setupFn);
        end
        
        %------------------------------------------------------------------
        
        function p = cleanup(p, str)
            if nargin < 2
                str = '';
            end
            if strcmp(str, 'force')
                rmdir(p.outputPath, 's');
            else
                if numel(dir(p.outputPath)) > 2
                    h = questdlg(['Really delete non-empty folder: ', p.outputPath, '?'], ...
                        'Deleting folder', 'Yes', 'No', 'Abort', 'No');
                    switch h
                        case 'Yes'
                            rmdir(p.outputPath, 's');
                        case 'Abort'
                            return
                    end
                elseif isfolder(p.outputPath)
                    rmdir(p.outputPath);
                end
            end
            p.controlList       = {};
            p.evaluationFolders = {};
            delete(p);
        end
        
        %------------------------------------------------------------------
        
        function [meanvals, vals] = loadAllFiles(p, fnm, field, warn)
            if nargin < 4
                warn = true;
            end
            emptyOutput = true;
            [nr, nc] = deal(numel(p.realizations), numel(p.controlList));
            meanvals = repmat({nan}, 1, nc);
            vals     = repmat({nan}, nr, nc);
            for kc = 1:nc
                nm = fullfile(p.evalOutputPath(kc), fnm);
                if isfile(nm)
                    meanvals{kc} = load(nm);
                    if ~isempty(field)
                        meanvals{kc} = meanvals{kc}.(field);
                    end
                end
                for kr = 1:nr
                    realNo = p.realizations(kr);
                    nm = fullfile(p.memberOutputPath(kc,realNo), fnm);
                    if isfile(nm)
                        vals{kr,kc} = load(nm);
                        emptyOutput = false;
                        if ~isempty(field)
                            vals{kr,kc} = vals{kr,kc}.(field);
                        end
                    end
                end
            end
            if emptyOutput && warn
                warning('No files: %s was currently found', fnm);
            end
        end
        
        %------------------------------------------------------------------
        
        function control = schedule2control(p, schedule)
            % map schedule to scaled control-vector
            scale = @(v, bnds)(v-bnds(1))/(bnds(2)-bnds(1));
            control = zeros(numel(p.maps.type), 1);
            for k = 1:numel(control)
                [sno, wno, bnds, tp] = deal(p.maps.stepNo(k), p.maps.wellNo(k), p.maps.bounds(k,:), p.maps.type{k});
                if p.maps.isTarget(k)
                    assert(strcmp(tp, schedule.control(sno).W(wno).type), 'Schedule contains non-expected well type');
                    control(k) = scale(schedule.control(sno).W(wno).val, bnds);
                else
                    control(k) = scale(schedule.control(sno).W(wno).lims.(tp), bnds);
                end
            end
        end
        
        %------------------------------------------------------------------
        
        function schedule = control2schedule(p, control)
            % map scaled control-vector to schedule
            schedule = p.initialGuess;
            scaleBack = @(val, bnds) bnds(1) + val*(bnds(2)-bnds(1));
            for k = 1:numel(control)
                [sno, wno, bnds, tp] = deal(p.maps.stepNo(k), p.maps.wellNo(k), p.maps.bounds(k,:), p.maps.type{k});
                if p.maps.isTarget(k)
                    assert(strcmp(tp, schedule.control(sno).W(wno).type), 'Schedule contains non-expected well type');
                    schedule.control(sno).W(wno).val = scaleBack(control(k), bnds);
                    % also update if part of lims
                    if isfield(schedule.control(sno).W(wno).lims, tp)
                        schedule.control(sno).W(wno).lims.(tp) = scaleBack(control(k), bnds);
                    end
                else
                    schedule.control(sno).W(wno).lims.(tp) = scaleBack(control(k), bnds);
                end
            end
        end
        
        %------------------------------------------------------------------
        
        function opt = getObjectiveOptions(p, obj)
            if isequal(obj, p.objective)
                opt = p.objectiveOpts;
            else
                opt = {};
            end
        end
        
        % -----------------------------------------------------------------
        
        function p = refreshControlList(p)
            d     = dir(p.outputPath);
            fldrs = regexp({d.name}, '\d+', 'match');
            p.evaluationFolders = [fldrs{:}];
            % load previously computed controls
            p.controlList = cell(0);
            for k = 1:numel(p.evaluationFolders)
                controlFile = fullfile(p.outputPath, p.evaluationFolders{k}, 'control.mat');
                if isfile(controlFile)
                    tmp = load(controlFile);
                    p.controlList{k} = tmp.control;
                end
            end
        end
    end
end

% Other functions
function f = checkFunction(f)
if ~isa(f, 'function_handle')
    assert(isa(f, 'char'), 'Expected function or string, but got %s', class(f));
    f = str2func(f);
end
end

function g = mapGradient(grad, maps)
scale = @(gi, bnds)gi*(bnds(2)-bnds(1));
g = zeros(numel(maps.type), 1);
for k = 1:numel(g)
    [sno, wno, bnds, tp] = deal(maps.stepNo(k), maps.wellNo(k), maps.bounds(k,:), maps.type{k});
    g(k) = scale(grad.(tp)(wno, sno) ,bnds);
end
end

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
