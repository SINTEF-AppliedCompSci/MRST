classdef ModelEnsemble < handle
    properties 
        nMembers    =  nan;
        name   
        setupFn     =  '';
        diagnostics =  [];  
        simulations =  [];
        evaluations =  []; % for optimization purposes
        models      =  [];
        writeToDisk   = true;
        writeModelsToDisk = true;
        singleScenario    = true;  
        plotProgress  = true;
        limitedOutput = false;
        nWorkers = 0;
        runBackground = false;
        figures
    end
    properties
        outputPath
    end
    properties (SetAccess = private)
        controlList    = [];
        outNms         = struct('sim',   'simulation', ...
                                'diagn', 'diagnostics', ...
                                'obj',   'objective');
    end
    
    methods
        function ensemble = ModelEnsemble(name, varargin)
            % ensemble = ModelEnsemble(name, 'pn1', 'pv1', ...)
            % Create or load model ensemble
            % REQUIRED PARAMETERS
            % name : name of ensemble
            %
            % The following optional parameter is REQUIRED on creation of ensemble
            % setupFn  - function handle to setup-function. Should take as
            %            input realization number and should output struct
            %            with fields 'model', 'W', and optionally 'state0',
            %            'schedule' and 'simOpts'
            %
            % OPTIONAL PARAMETERS
            % nMembers       - number of realizations in ensemble (if not 
            %                  given, it will be set based on output of
            %                  setupFn
            % singleScenario - if true, only a single control and/or
            %                  schedule is considered (default true)
            % nWorkers       - number of external matlab sessions used for 
            %                  computations (default 0)
            % runBackground  - If true, all computations will be run in
            %                  background, if false and nWorkers > 0,
            %                  computation progress will be displayed
             
            [ensemble, extra] = merge_options(ensemble, varargin{:});
            ensemble.outputPath = fullfile(mrstOutputDirectory(), name);
            
            % cleanup - option
            opt = struct('cleanup', false);
            opt = merge_options(opt, extra{:});

            if opt.cleanup
                cleanupOutputDirectory(ensemble);
            end
            
            % get/load/write setup
            ensemble.name = name;
            ensemble = handleSetup(ensemble);
                                 
            % setup handlers
            ensemble.setupSimulationHandlers();
            ensemble.setupDiagnosticsHandlers();
            
            [~, nm] = fileparts(ensemble.outputPath);
            ensemble.name = nm;
            
            % check if setupFn is of 'old' type
            ensemble.setupFn = checkSetupFn(ensemble.setupFn);
            
            %ensemble.setupModelHandlers();
%             NEEDS rewriting in terms of handlers
%             ensemble.setupEvaluationHandlers();
%             
%             if ~isempty(ensemble.evaluations)
%                 ensemble.controlList = cell(0);
%                 for k = 1:numel(ensemble.objFolders)
%                     controlFile = fullfile(pth, ensemble.objFolders{k}, 'controls.mat');
%                     if isfile(controlFile)
%                         tmp = load(controlFile);
%                         ensemble.controlList{k} = tmp.controls;
%                     end
%                 end
%             end
        end
        
        %------------------------------------------------------------------
        function setupSimulationHandlers(ensemble, simNo)
            if nargin == 1 % setup-stage search for existing folders
                d    = dir(ensemble.outputPath);
                fnms = {d.name};
                simfldrs = regexp(fnms, [ensemble.outNms.sim,'\d+'], 'match');
                simNo = 1:numel([simfldrs{:}]);
            end

            fields = {'ws', 'states', 'reports'};
            
            for kd = 1:numel(fields)
                for ks = 1:numel(simNo)
                    curSimNo = simNo(ks);
                    ddir = fullfile(ensemble.outputPath,  sprintf('simulation%4.4d', curSimNo));
                    if ~isfolder(ddir)
                        mkdir(ddir);
                    end
                    ensemble.simulations(ks).(fields{kd}) = arrayfun(@(memNo)...
                        ResultHandler('dataDirectory', ddir, 'dataFolder', sprintf('realization%4.4d', memNo), ...
                        'writeToDisk', ensemble.writeToDisk, 'dataPrefix', fields{kd}), ...
                        (1:ensemble.nMembers), 'UniformOutput', false);
                end
            end
        end
        
        %------------------------------------------------------------------
        function setupDiagnosticsHandlers(ensemble, diagnNo, fields)
            if nargin == 1 % setup-stage search for existing folders
                d    = dir(ensemble.outputPath);
                fnms = {d.name};
                diagnfldrs = regexp(fnms, [ensemble.outNms.diagn,'\d+'], 'match');
                diagnfldrs = [diagnfldrs{:}];
                diagnNo = 1:numel(diagnfldrs);
                
                if isempty(diagnNo)
                    return;
                else
                   fields = {'states', 'D', 'WP', 'RTD', 'other','statistics'};
                end
            elseif nargin == 2
                fields = {'states', 'D', 'WP', 'RTD', 'other','statistics'};
            end
            
            if ~isempty(fields)
                for kd = 1:numel(diagnNo)
                    curDiagnNo = diagnNo(kd);
                    for kf = 1:numel(fields)
                        ensemble.diagnostics(curDiagnNo).(fields{kf}) = ...
                            ResultHandler('dataDirectory', ensemble.outputPath, ...
                            'dataFolder', sprintf('diagnostics%4.4d', curDiagnNo), ...
                            'writeToDisk', ensemble.writeToDisk, ...
                            'dataPrefix', fields{kf} );
                    end
                end
            end
        end
        
        %------------------------------------------------------------------
        
        function str = simOutputPath(ensemble, simNo)
            str = fullfile(ensemble.outputPath, ...
                           sprintf('simulation%4.4d', simNo));
        end
        
        %------------------------------------------------------------------
        
        function str = diagnOutputPath(ensemble, diagnNo)
            str = fullfile(ensemble.outputPath, ...
                           sprintf('diagnostics%4.4d', diagnNo));
        end
        
        %------------------------------------------------------------------
        
        function str = memberOutputPath(ensemble, simNo, memNo)
            str = fullfile(ensemble.simOutputPath(simNo), ...
                           sprintf('realization%4.4d', memNo) );
        end
        
        %------------------------------------------------------------------
            
        function computeDiagnostics(ensemble, varargin)
            opt = struct('outputDir',                              [], ...
                         'control',                                [], ...
                         'memberIx',            (1:ensemble.nMembers), ...
                         'rtd',                            'compute',  ...
                         'diagnNo',                                [], ...
                         'nWorkers',                ensemble.nWorkers, ...
                         'runBackground',      ensemble.runBackground);
                     
            [opt, extra] = merge_options(opt, varargin{:});
            
            if ensemble.singleScenario
                diagnNo = 1;
            elseif ~isempty(opt.diagnNo)
                diagnNo = opt.diagnNo;
            else
                diagnNo = numel(ensemble.diagnostics)+1;
            end
            
            % check if some or all already have been computed
            stat = ensemble.getDiagnosticsStatus(diagnNo);
            if ~isempty(stat.fields)
                isComputed = stat.status(:, strcmp(stat.fields, 'RTD'));
                opt.memberIx = opt.memberIx(~isComputed(opt.memberIx));
                if isempty(opt.memberIx)
                    fprintf('Diagnostics already computed\n');
                    return;
                end
            end
            
            if isempty(opt.outputDir)
                outputDir = ensemble.diagnOutputPath(diagnNo);
            end
            
            ensemble.setupDiagnosticsHandlers(diagnNo);            
            
            % check if control-file already exists for scenario
            cntNm = fullfile(outputDir, 'control.mat');
            if ~isfile(cntNm)
                dispif(mrstVerbose, 'Initial setup of diagnostics computation for scenario %d ...', diagnNo);
                tmp = ensemble.setupFn(1);
                control = struct('W', tmp.W);
                save(cntNm, 'control');
                dispif(mrstVerbose, 'done.')
            else
                tmp = load(cntNm);
                control = tmp.control;
            end
            
            [nWork, runBack] = deal(opt.nWorkers, opt.runBackground);
            nWork = max(runBack, nWork); % if background, use at least one worker
            
            inputs = [{'outputDir', outputDir, 'rtd', 'compute'}, extra];
            if nWork == 0
                inputs = [inputs, {'showWaitbar', false}];
                computeEnsembleDiagnostics(ensemble.setupFn, 'memberIx', opt.memberIx, ...
                                           'control', control, inputs{:}); 
            else % launch new matlab sessions
                inputs = [inputs, {'showWaitbar', false}];
                argFuncs = struct('name',              'memberIx', ...
                                  'n',        numel(opt.memberIx), ...
                                  'fn',        @(i)opt.memberIx(i));
                if ~runBack % run sessions with with progress plot
                    ensemble.launchSessions(@computeEnsembleDiagnostics, ...
                        [{ensemble.setupFn, 'control', cntNm}, inputs], argFuncs, ...
                        outputDir, nWork);
                else %run all sessions in background
                    ensemble.launchSessionsBackground(@computeEnsembleDiagnostics, ...
                        [{ensemble.setupFn, 'control', cntNm}, inputs], argFuncs, ...
                        outputDir, nWork);
                end
            end
        end
        
        %------------------------------------------------------------------
            
        function runSimulations(ensemble, varargin)
            opt = struct('outputDir',                              [], ...
                         'schedule',                               [], ...
                         'memberIx',            (1:ensemble.nMembers), ...
                         'simNo',                                  [], ...
                         'nWorkers',                ensemble.nWorkers, ...
                         'runBackground',      ensemble.runBackground);
                     
            [opt, extra] = merge_options(opt, varargin{:});
            
            if ensemble.singleScenario
                simNo = 1;
            elseif ~isempty(opt.diagnNo)
                simNo = opt.diagnNo;
            else
                simNo = numel(ensemble.diagnostics)+1;
            end
            
            if isempty(opt.outputDir)
                outputDir = ensemble.simOutputPath(simNo);
            end
            
            ensemble.setupSimulationHandlers(simNo)            
    
            % check if schedule-file already exists for scenario
            schNm = fullfile(outputDir, 'schedule.mat');
            if ~isfile(schNm)
                dispif(mrstVerbose, 'Initial setup of simulations for scenario %d ...', simNo);
                tmp = ensemble.setupFn(1);
                schedule = tmp.schedule;
                %control = struct('W', tmp.W);
                save(schNm, 'schedule');
                dispif(mrstVerbose, 'done.')
            else
                tmp = load(schNm);
                schedule = tmp.schedule;
            end
            
            % check if some or all already have been computed
            stat = ensemble.getSimulationStatus(simNo);
            if ~isempty(stat.fields)
                nsteps = numel(schedule.step.val);
                isComputed = (stat.status(:, strcmp(stat.fields, 'ws')) == nsteps);
                opt.memberIx = opt.memberIx(~isComputed(opt.memberIx));
                if isempty(opt.memberIx)
                    fprintf('Simulations already run\n');
                    return;
                end
            end
            
            [nWork, runBack] = deal(opt.nWorkers, opt.runBackground);
            nWork = max(runBack, nWork); % if background, use at least one worker
            
            inputs = [{'outputDir', outputDir}, extra];
            if nargout > 0 || nWork == 0
                simulateEnsemble(ensemble.setupFn, 'memberIx', opt.memberIx, 'schedule', schedule, inputs{:});
            else % launch new matlab sessions
                argFuncs = struct('name',              'memberIx', ...
                                  'n',        numel(opt.memberIx), ...
                                  'fn',        @(i)opt.memberIx(i));
                if ~runBack % run sessions with with progress plot
                    ensemble.launchSessions(@simulateEnsemble, ...
                        [{ensemble.setupFn, 'schedule', schNm}, inputs], argFuncs, ...
                        outputDir, nWork);
                else %run all sessions in background
                    ensemble.launchSessionsBackground(@simulateEnsemble, ...
                        [{ensemble.setupFn, 'schedule', schNm}, inputs], argFuncs, ...
                        outputDir, nWork);
                end
            end    
        end
        
        %------------------------------------------------------------------
        
        function value = getMemberObjective(ensemble, objective, memberNo, varargin)
            [control, controlNo, curpath] = ensemble.handleControlInfo(memberNo, varargin{:});
            fileNm = fullfile(curpath, [func2str(objective), '.mat']);
            if isfile(fileNm) 
                tmp = load(fileNm);
                value = tmp.value;
            else
                if ~isfolder(curpath)
                    mkdir(curpath);
                end
                value = ensemble.evalFn(control, ensemble.setupFn, ...
                                'outputPath', ensemble.memberOutputPath(controlNo, memberNo), ...
                                'memberNo',   memberNo, ...
                                'objective',  objective); 
                save(fileNm, 'value');
            end
            if ensemble.plotProgress
                ax = getValidAxes(ensemble, 'main');
                hold(ax, 'on');
                if numel(value) == 1
                    plot(ax, controlNo, value, 'ob');
                    drawnow
                end
            end
        end
        
        %------------------------------------------------------------------
        
        function value = getObjective(ensemble, objective, varargin)
            [control, controlNo, curpath] = ensemble.handleControlInfo([], varargin{:});
            fileNm = fullfile(curpath, [func2str(objective), '.mat']);
            if isfile(fileNm)
                tmp = load(fileNm);
                value = tmp.value;
            else
                if ensemble.parallel > 0
                    % handle computations in seperate matlab sessions
                    ensemble.launchMultipleSessions(objective, 'control', control, 'controlNo', controlNo);
                end
                % compute or read individual values
                [value, numOK] = deal(0);
                for memberNo = 1:ensemble.nMembers
                    vi = ensemble.getMemberObjective(objective, memberNo, 'control', control, 'controlNo', controlNo);
                    if ~isempty(vi)
                        numOK = numOK +1;
                        value = value + vi;
                    end
                end
                value = value/numOK;
                save(fileNm, 'value');
            end
            if ensemble.plotProgress
                ax = getValidAxes(ensemble, 'main');
                hold(ax, 'on');
                if numel(value) == 1
                    plot(ax, controlNo, value, 'xr');
                    drawnow
                end
            end
        end
        
        %------------------------------------------------------------------
        
        function gradient = getMemberGradient(ensemble, objective, memberNo, varargin)
            [control, controlNo, curpath] = ensemble.handleControlInfo(memberNo, varargin{:});
            fileNm = fullfile(curpath, ['gradient_', func2str(objective), '.mat']);
            if isfile(fileNm)
                tmp = load(fileNm);
                gradient = tmp.gradient;
            else
                if ~isfolder(curpath)
                    mkdir(curpath);
                end
                [~, gradient] = emsemble.evalFn(control, ensemble.setupFn, ...
                                    'outputPath', ensemble.memberOutputPath(controlNo, memberNo), ...
                                    'memberNo',   memberNo, ...
                                    'objective',  objective);
                save(fileNm, 'gradient');
            end
        end
        
        %------------------------------------------------------------------
        
        function gradient = getGradient(ensemble, objective, varargin)
            [control, controlNo, curpath] = ensemble.handleControlInfo([], varargin{:});
            fileNm = fullfile(curpath, ['gradient_', func2str(objective), '.mat']);
            if isfile(fileNm)
                tmp = load(fileNm);
                gradient = tmp.gradient;
            else
                if esemble.parallel > 0
                    % handle computations in seperate matlab sessions
                    ensemble.launchMultipleSessions(objective, 'control', control, 'controlNo', controlNo);
                end
                % compute or read individual values
                [gradient, numOK] = deal([], 0);
                for memberNo = 1:ensemble.nMembers
                    gi = ensemble.getMemberGradient(objective, memberNo, 'control', control, 'controlNo', controlNo);
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
        
        function values = getComputedObjectives(ensemble, objective)
            % loop through schedule-list and collect all objective values
            fileNm = [func2str(objective), '.mat'];
            ns = numel(ensemble.scheduleList);
            ne = ensemble.nMembers;
            values = struct('total', nan(1, ns), 'member', nan(ne, ns));
            for ks = 1:ns
                [~, curpath, ~] = ensemble.proccessScheduleInfo([], 'controlNo', ks);
                outNm = fullfile(curpath, fileNm);
                if isfile(outNm)
                    tmp = load(outNm);
                    values.total(ks) = tmp.value;
                end
                for ke = 1:ne
                    [~, curpath, ~] = ensemble.proccessScheduleInfo(ke, 'controlNo', ks);
                    outNm = fullfile(curpath, fileNm);
                    if isfile(outNm)
                        tmp = load(outNm);
                        values.member(ke, ks) = tmp.value;
                    end
                    
                end
            end         
        end
        
        %------------------------------------------------------------------
    
        function launchSessions(ensemble, fn, args, argFuncs, outDir, nWork)
            if nargin < 6
                nWork = ensemble.nWorkers;
            end
            nComput = argFuncs(1).n;
            [startTimes, endTimes] = deal(nan(nComput, 1));
            activePrev = false(nComput, 1);

            [compNo, numActive] = deal(0);
            while any(isnan(endTimes))
                [numActive, activeIx] = getNumberOfSessions(outDir, nComput);
                isDone = and(activePrev, ~activeIx);
                activePrev = activeIx;
                endTimes(isDone) = datenum(clock);
                
                if numActive < nWork && compNo < nComput
                    % launch new
                    compNo  = compNo + 1;
                    %numActive = numActive + 1;
                    startTimes(compNo) =  datenum(clock);
                    curargs = args;
                    for ka = 1:numel(argFuncs)
                        curargs = [curargs, {argFuncs(ka).name, argFuncs(ka).fn(compNo)}]; %#ok
                    end
                    progFn = fullfile(outDir, sprintf('progress%d.mat', compNo));
                    fprintf('Starting new Matlab session ...\n')
                    evalFunWrapper(fn, curargs, 'progressFileNm', progFn, 'exitWhenDone', true);
                else
                    pause(1);
                end
                if ensemble.plotProgress
                    plotSessionProgress(ensemble, startTimes, endTimes)
                end
            end
        end
        
        %------------------------------------------------------------------
        
        function launchSessionsBackground(ensemble, fn, args, argFuncs, outDir, nWork)
            if nargin < 6
                nWork = max(1, ensemble.nWorkers);
            end
            nComput = argFuncs(1).n;
            binSize = ceil(nComput/nWork);
            for k = 1:nWork
                bin = (1+(k-1)*binSize) : min(k*binSize,nComput);
                if ~isempty(bin)
                    curargs = args;
                    for ka = 1:numel(argFuncs)
                        curargs = [curargs, {argFuncs(ka).name, argFuncs(ka).fn(bin)}]; %#ok
                    end
                    fprintf('Starting new Matlab session ...\n')
                    evalFunWrapper(fn, curargs, 'exitWhenDone', true);
                end
            end
        end
        
        %------------------------------------------------------------------
        
%         function launchMultipleSessions(ensemble, objective, varargin)
%             [~, controlNo, curpath] = ensemble.handleControlInfo([], varargin{:});
%             controlFile = fullfile(curpath, 'control.mat');
%             assert(isfile(controlFile), 'Unable to locate control-file');
%             memberPathList = cell(ensemble.nMembers, 1);
%             for memberNo = 1:ensemble.nMembers
%                 memberPathList{memberNo} = ensemble.memberOutputPath(controlNo, memberNo);
%             end
%             if ensemble.plotProgress
%                 %ensemble.figures.sessions = figure;
%                 ax = ensemble.getValidAxes('sessions');
%                 plot(ax, datetime(now*ones(2, ensemble.nMembers), 'convertFrom', 'datenum'), ...
%                          repmat(1:ensemble.nMembers, [2 ,1]), 'LineWidth', 10, 'Color', 'g');
%                 ax.FontSize = 12;
%                 ax.XLim  = datetime([now-.01, now], 'convertFrom', 'datenum');
%                 title(ax, ['Computing ', func2str(objective) ' for control number ', num2str(controlNo) ]);
%             end
%             
%             [memberNo, numActive] = deal(0);
%             [startTimes, endTimes] = deal(nan(ensemble.nMembers, 1));
%             activePrev = false(ensemble.nMembers, 1);
%             while any(isnan(endTimes))%numActive < ensemble.parallel && memberNo < ensemble.nMembers
%                 [numActive, activeIx] = getNumberOfSessions(memberPathList);
%                 isDone = and(activePrev, ~activeIx);
%                 activePrev = activeIx;
%                 endTimes(isDone) = datenum(clock);
%                     
%                 if numActive < ensemble.parallel && memberNo < ensemble.nMembers 
%                     % launch new
%                     memberNo  = memberNo + 1;
%                     %numActive = numActive + 1;
%                     startTimes(memberNo) =  datenum(clock); 
%                     % locate control-file
%                     fileNm = fullfile(curpath, 'control.mat');
%                     assert(isfile(fileNm));
%                     if ~isfolder(memberPathList{memberNo})
%                         mkdir(memberPathList{memberNo});
%                     end
%                     evaluateObjectiveWrapper(ensemble.evalFn, fileNm, ...
%                         ensemble.setupFn, 'outputPath', memberPathList{memberNo}, ...
%                         'memberNo', memberNo, 'objective', objective);
%                     %system(str)
%                     % here we should check that first progress file is
%                     % written before continuing
%                 else
%                     pause(1);
%                 end
%                 if ensemble.plotProgress
%                     plotSessionProgress(ensemble, startTimes, endTimes)
%                 end
%             end
%             
%             %if ensemble.plotProgress
%             %    close(ensemble.figures.sessions);
%             %end
%         end
                
        %------------------------------------------------------------------
        
        function [control, controlNo, curpath] = handleControlInfo(ensemble, memberNo, varargin)
            opt = struct('control', [], 'controlNo', []);
            opt = merge_options(opt, varargin{:});
            [control, controlNo] = deal(opt.control, opt.controlNo);
            if isempty(controlNo) % we need to check exisiting controls
                assert(~isempty(control), 'Need either control or controlNo ...')
                % search for mathcing previous controls (use single precision)
                for k = numel(ensemble.controlList):-1:1
                    if norm( single(control) - single(ensemble.controlList{k})) == 0
                        controlNo = k;
                        break;
                    end
                end
            end
            
            % if no controlNo is found, we have a new control, record accordingly
            if isempty(controlNo) 
                controlNo   = numel(ensemble.controlList) +1;
                controlFolder = ensemble.simOutputPath(controlNo);
                if ~isfolder(controlFolder)
                    mkdir(controlFolder);
                end
                % save control and add to control-list
                save(fullfile(controlFolder, 'control.mat'), 'control');
                ensemble.controlList = [ensemble.controlList, {control}];
            end
            % if schedule is empty, attempt to load it
            controlFolder = ensemble.simOutputPath(controlNo);
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
                curpath = ensemble.memberOutputPath(controlNo, memberNo);
            end
        end
        
        %------------------------------------------------------------------
        
        function ax = getValidAxes(ensemble, figname)
            if ~isfield(ensemble.figures, figname) || ~isvalid(ensemble.figures.(figname))
                ensemble.figures.(figname) = figure('ToolBar', 'none', 'MenuBar', 'none');
            end
            fig = ensemble.figures.(figname);
            % search for axes
            isAx = arrayfun(@(x)isa(x, 'matlab.graphics.axis.Axes'), fig.Children);
            if any(isAx)
                ax = fig.Children(find(isAx, 1, 'last'));
            else
                ax = axes(fig);
            end
        end
        
        %------------------------------------------------------------------
        
        function v = getDiagnosticsStatus(ensemble, diagnNo)
            if nargin < 2 && ensemble.singleScenario
                diagnNo = 1;
            end
            v = struct('fields', {{}}, 'status', false);
            if ~isempty(ensemble.diagnostics) && numel(ensemble.diagnostics) >= diagnNo
                diagn = ensemble.diagnostics(diagnNo);
                flds = fieldnames(diagn);
                status = false(ensemble.nMembers, numel(flds));
                for k = 1:numel(flds)
                    status(diagn.(flds{k}).getValidIds, k) = true;
                end
                [v.fields, v.status] = deal(flds(:)', status);
            end
                   
        end
        
        %------------------------------------------------------------------
        
        function v = getSimulationStatus(ensemble, simNo)
            if nargin < 2 && ensemble.singleScenario
                simNo = 1;
            end
            v = struct('fields', {{}}, 'status', false);
            if ~isempty(ensemble.simulations) && numel(ensemble.simulations) >= simNo
                sim = ensemble.simulations(simNo);
                flds = fieldnames(sim);
                status = zeros(ensemble.nMembers, numel(flds));
                for k = 1:numel(flds)
                    nr = numel(sim.(flds{k}));
                    status(1:nr, k) = cellfun(@(x)numel(x.getValidIds), sim.(flds{k}));   
                end
                [v.fields, v.status] = deal(flds(:)', status);
            end
        end
        
        %------------------------------------------------------------------
        
        function viewStatus(ensemble, scenarioNum)
            if nargin < 2 && ensemble.singleScenario
                scenarioNum = 1;
            end
            [data, names] = deal({});
            ds = ensemble.getDiagnosticsStatus(scenarioNum);
            if any(strcmp(ds.fields, 'RTD')) 
                data  = [data, {ds.status(:, strcmp(ds.fields, 'RTD'))}];
                names = [names, {'diagnostics'}];
            end
            ss = ensemble.getSimulationStatus(scenarioNum);
            if any(strcmp(ss.fields, 'ws'))
                outputDir = ensemble.simOutputPath(scenarioNum);
                schNm = fullfile(outputDir, 'schedule.mat');
                tmp = load(schNm);
                nsteps = numel(tmp.schedule.step.val);
                data  = [data, {100*ss.status(:, strcmp(ss.fields, 'ws'))/nsteps}];
                names = [names, {'simulations'}];
            end
            
            plotStatus(data, names);
        end
    end
end

%--------------------------------------------------------------------------
% PRIVATE HELPERS
%--------------------------------------------------------------------------

function nm = getValidFunctionName(fn)
assert(isa(fn, 'function_handle'));
nm = func2str(fn);
if ~isvarname(nm)  
    % contains special characters, attempt to extract function name in
    % anonymous expression
    nm = regexp(nm, ')(\w*)(', 'once', 'tokens');
    while isa(nm, 'cell')
        nm = nm{1};
    end
end
end

%------------------------------------------------------------------
% function schedule = getScheduleForCurrentWell(schedule, W)
% ns = schedule.step.control(end);
% assert(ns==1 || ~any(diff(arrayfun(@(x)numel(x.W), schedule.control)) ~= 0) , ...
%        'Schedule appears to be non-uniform');
% assert(numel(W) == numel(schedule.control(1).W), 'Non-compatible schedules') 
% flds = {'type', 'val', 'lims'};
% for ks = 1:ns
%     tmp = schedule.control(ks).W;
%     schedule.control(ks).W = W;
%     curflds = flds(isfield(tmp, flds));
%     for kw = 1:numel(W)
%         for kf = 1:numel(curflds)
%             schedule.control(ks).W(kw).(curflds{kf}) = tmp(kw).(curflds{kf});
%         end
%     end
% end
% end

%--------------------------------------------------------------------------

% function s = getCallString(fun, fileNm, obj, setupFn, memNo, pth)
% fun     = func2str(fun);
% obj     = func2str(obj);
% setupFn = func2str(setupFn);
% memNo   = num2str(memNo);
% 
% dsp = @(str)['''', str, ''''];
% 
% inp = [dsp(fileNm), ',', dsp(setupFn),      ',', ...
%        dsp('memberNo'),          ',', memNo,     ',', ...
%        dsp('objective'),         ',', dsp(obj),  ',', ...
%        dsp('outputPath'),        ',', dsp(pth),  ',', ...
%        dsp('writeProgressFile'), ',', '1'];
% 
% funcall = [fun, '(', inp, ')'];
% 
% %s = ['matlab -nodisplay -nosplash -nodesktop -r ', ...
% %    '"', funcall, ';",&']; 
%  
% s = ['matlab -nodisplay -nosplash -nodesktop -r ', ...
%      '"', funcall, ';exit",&'];   
% end

%------------------------------------------------------------------

function plotSessionProgress(ensemble, t1, t2)
ax = ensemble.getValidAxes('sessions');
if isempty(ax.Children)
    ax.Parent.Color = [1 1 1];
    ax.Color = [1 1 1];
    plot(ax, datetime(now*ones(2, ensemble.nMembers), 'convertFrom', 'datenum'), ...
         repmat(1:ensemble.nMembers, [2 ,1]), 'LineWidth', 15, 'Color', 'b');
    ax.FontSize = 14;
    ax.XLim  = datetime([now-.01, now], 'convertFrom', 'datenum');
    ax.YLabel.String = 'Job No';
    title(ax, 'Computing status');
end

doUpdate = true;
if doUpdate
    nMem = numel(t1);
    t = [t1, t2];
    t(and(~isnan(t1),isnan(t2)), 2) = now;
    for k = 1:nMem
        ax.Children(ensemble.nMembers+1-k).XData = datetime(t(k,:), 'convertFrom', 'datenum');
        ax.Children(ensemble.nMembers+1-k).YData = [k k];
    end
    ax.YTick = 1:nMem;
    axis(ax, 'tight')
    ax.YLim  = [.5, nMem+.5];
    ax.FontSize = 12;
    ax.XTickLabelRotation = 30;
end
end

%--------------------------------------------------------------------------

function plotStatus(data, titles)
nd = numel(data);
if nargin < 2
    titles = repmat({''}, [1, nd]);
end
if nd > 0
    f = limitedToolbarFigure;
    for k = 1:nd
        cur = data{k};
        ax = subplot(1, nd, k);
        if islogical(cur)
            ix = find(cur);
            plot(ax, ones(numel(ix),1), ix, 'LineStyle', 'none', 'Marker', 's', 'MarkerSize', 5, 'LineWidth', 4, 'MarkerEdgeColor', [0 .5 0]);
            ax.XAxis.Visible = 'off';
            axis equal tight
            ax.YLim = [0 numel(cur)];
        else
            plot(ax, [zeros(1, numel(cur)); cur(:)'], repmat((1:numel(cur)), [2, 1]), 'LineStyle', '-', 'LineWidth', 4, 'Color', [0 .5 0]);
            ax.XLim = [0 100];
        end
        ax.Title.String = titles{k};
        ax.YTick = 1:numel(cur);
    end
end
end       

%------------------------------------------------------------------

function [na, ix] = getNumberOfSessions(outdir, n)
ix = false(n,1);
fn = dir(fullfile(outdir, 'progress*.mat'));
na = numel(fn);
nms = {fn.name};

ii = regexp(nms, '\d+', 'match');
if ~isempty(ii)
    ii = cellfun(@str2num, [ii{:}]);
else
    ii = [];
end
ix(ii) = true;
end

%------------------------------------------------------------------

function  ensemble = handleSetup(ensemble)
setupFile  = fullfile(ensemble.outputPath, 'setup.mat');
if isempty(ensemble.setupFn)% previous setup must exist
    if ~(isfolder(ensemble.outputPath) && isfile(setupFile))
        error('Unable to find setup-details for ensemble: %s', ensemble.name);
    end
    tmp = load(fullfile(ensemble.outputPath, 'setup.mat'));
    ensemble.nMembers = tmp.nMembers;
    ensemble.setupFn  = tmp.setupFn;
end

if isempty(ensemble.nMembers)
    error('Number of ensemble members unset');
end

if ischar(ensemble.setupFn)
    assert(isfile(ensemble.setupFn), 'Unable to locate: %s', ensemble.setupFn);
    [pth, nm] = fileparts(ensemble.setupFn);
    if isempty(which(nm))
        fprintf('Setup function (%s.m) not on current path.\nAdding %s to path.\n', nm, pth);
        addpath(pth);
    end
    ensemble.setupFn = str2func(nm);
end

assert(isa(ensemble.setupFn, 'function_handle'), 'Setup function of unreckognized type');
% don't support anonymous functions
f = functions(ensemble.setupFn);
assert(~strcmp(f.type, 'anonymous'), 'Setup function of anonymous type is not supported');
assert(~isempty(f.file), 'Setup file location not found');

% write setup info to outputPath
if ~isfolder(ensemble.outputPath)
    mkdir(ensemble.outputPath);
end

nMembers = ensemble.nMembers;                   %#ok
setupFn  = f.file;                              %#ok
save(setupFile, 'nMembers', 'setupFn')
end

%------------------------------------------------------------------

function cleanupOutputDirectory(ensemble)
if numel(dir(ensemble.outputPath)) > 2
    h = questdlg(['Really delete non-empty folder: ', ensemble.outputPath, '?'], ...
        'Deleting folder', 'Yes', 'No', 'Abort', 'No');
    switch h
        case 'Yes'
            rmdir(ensemble.outputPath, 's');
        case 'Abort'
            return
    end
elseif isfolder(ensemble.outputPath)
    rmdir(ensemble.outputPath);
end
end

%--------------------------------------------------------------------------
% backward compatibility, remove at some time -----------------------------
function f = checkSetupFn(f)
if nargout(f) > 1
    warning('Format of setupFn should be updated')
    f = @(n)wrapSetupFn(f, n);
end
end

function out = wrapSetupFn(f, n)
[state0, schedule] = deal([]);
if nargout(f) == 2
    [model, W] = f(n);
else
    [model, W, state0] = f(n);
end
out = struct('model', model, 'W', W, 'state0', state0, 'schedule', schedule);
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
