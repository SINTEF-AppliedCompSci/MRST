classdef MRSTRegressionTest

    properties
        group = 'general'
        name
        setup
        setupArgs = {}
        problemInput = {}
        compareStates = true
        compareWellSols = true
        compareReports = true
        tolerance = struct('wellSols', 0, ...
                           'states'  , 0, ...
                           'reports' , 0)
        verbose = true
        compareFn = []
    end

    methods
        function rt = MRSTRegressionTest(setup, varargin)
            if nargin == 0
                return
            end

            rt.setup = setup;
            [rt, ~] = merge_options(rt, varargin{:});

            if isempty(rt.name)
                rt.name = defaultName(setup);
            end
            if isempty(rt.compareFn)
                rt.compareFn = @comparePackedRegressionResults;
            end
        end

        function report = runRegressionTest(rt, varargin)
            opt = struct('deletePrevious', true);
            opt = merge_options(opt, varargin{:});

            setup = rt.resolveSetup();
            directory = rt.getDataDirectory();
            if ~isfolder(directory)
                mkdir(directory);
            end

            hash = rt.getRegressionTestHash(setup);
            time = rt.getTimeString();
            namePrev = rt.findExisting(hash, time);

            if ~isempty(namePrev)
                problemPrev = rt.buildProblem(setup, namePrev);
            else
                warning(['No existing results for test case %s. Results ' , ...
                         'of this test will be inconclusive.'], rt.name);
                problemPrev = [];
            end

            nameCurr = [hash, '_', time];
            problemCurr = rt.buildProblem(setup, nameCurr);

            if rt.verbose
                rt.printHeader();
            end

            try
                simulatePackedProblem(problemCurr);
                details = [];
            catch ex
                details = ex.message;
            end

            if ~isempty(namePrev)
                report = rt.compareFn(rt, problemCurr, problemPrev, details);
                if report.passed == 1 && opt.deletePrevious
                    rmdir(fullfile(directory, namePrev), 's');
                end
            else
                if ~isempty(details), details = [details, '. ']; end
                details = [details, 'No existing results to compare against'];
                report = struct('passed', -1, 'details', details);
            end

            if rt.verbose
                rt.printFooter(report.passed);
            end

            save(fullfile(directory, ['regreport_', nameCurr, '.mat']), 'report');
        end

        function reports = getTestReports(rt)
            directory = rt.getDataDirectory();
            if ~isfolder(directory)
                reports = {};
                return
            end
            listing = dir(directory);
            reports = cell(numel(listing), 1);
            for i = 1:numel(listing)
                if listing(i).isdir || ~contains(listing(i).name, 'regreport_')
                    continue
                end
                data = load(fullfile(directory, listing(i).name));
                reports{i} = data.report;
            end
            reports = reports(~cellfun(@isempty, reports));
        end
    end

    methods (Access = private)
        function setup = resolveSetup(rt)
            if isa(rt.setup, 'function_handle')
                setup = rt.setup(rt.setupArgs{:});
            elseif ischar(rt.setup) || isstring(rt.setup)
                setup = feval(char(rt.setup), rt.setupArgs{:});
            else
                setup = rt.setup;
            end
        end

        function problem = buildProblem(rt, setup, name)
            if isPackedProblem(setup)
                problem = copyPackedProblem(setup, rt.group, ...
                    'Directory', rt.getDataDirectory(), ...
                    'Name', name);
                if ~isempty(rt.problemInput)
                    warning('mrst:regressionTestProblemInputIgnored', ...
                        'problemInput is ignored for packed-problem setups.');
                end
                return
            end

            packArgs = {'Directory', rt.getDataDirectory(), ...
                        'Name', name};
            if isfield(setup, 'description') && ~isempty(setup.description)
                packArgs = [packArgs, {'Description', setup.description}];
            end
            if isfield(setup, 'modules') && ~isempty(setup.modules)
                packArgs = [packArgs, {'Modules', setup.modules}];
            end
            problem = packSimulationProblem(setup.state0, setup.model, setup.schedule, ...
                rt.group, packArgs{:}, rt.problemInput{:});
        end

        function hash = getRegressionTestHash(rt, setup)
            if isPackedProblem(setup)
                setup = setup.SimulatorSetup;
            end
            hash = obj2hash({setup.state0, ...
                             obj2hash(setup.model, 'skip', 'inputdata'), ...
                             setup.schedule});
            if isempty(rt.problemInput)
                return
            end
            piHash = cellfun(@(ip) obj2hash(ip), rt.problemInput, 'UniformOutput', false);
            hash = obj2hash(strjoin([hash, piHash], '_'));
        end

        function directory = getDataDirectory(rt)
            directory = fullfile(mrstOutputDirectory(), 'reg-tests', rt.group, rt.name);
        end

        function existing = findExisting(rt, hash, timeCurr)
            directory = rt.getDataDirectory();

            existing = [];
            if ~exist(directory, 'dir'); return; end
            dirNames = ls(directory);
            if isempty(dirNames), existing = []; return; end
            if size(dirNames, 1) > 1
                pad = repmat(' ', size(dirNames, 1), 1);
                dirNames = reshape([dirNames, pad]', 1, []);
            end
            fmt = rt.dateTimeFormat();

            dirNames = regexp(char(dirNames), '\s+', 'split');
            nd       = numel(dirNames);
            existing = cell(nd, 1);
            for i = 1:nd
                if ~contains(dirNames{i}, hash) || ...
                   ~isfolder(fullfile(directory, dirNames{i}))
                    continue;
                end
                existing{i} = dirNames{i};
            end
            existing = existing(~cellfun(@isempty, existing));
            if isempty(existing), return; end
            timeCurr = datetime(timeCurr, 'Format', fmt);
            timePrev = cellfun(@(ext) datetime(ext(end-numel(fmt)+1:end), ...
                                               'Format', fmt), existing);
            if numel(timePrev) > 1
                warning(['More than one set of existing results found. ' , ...
                         'Using most recent results for regression testing.']);

                [timePrev, ix] = sort(timePrev);
                existing       = existing(ix);

                prompt = sprintf('Do you want to delete older results? y/n [y]: ');
                str    = input(prompt, 's');
                if strcmpi(str, 'n')
                    fprintf('Ok, will not remove results.\n');
                else
                    for i = 1:numel(existing)-1
                        fn = fullfile(directory, existing{i});
                        rmdir(fn, 's');
                    end
                end
            end
            timePrev = timePrev(end);
            assert(timePrev(end) < timeCurr, ...
                ['Previous results are stored with timestamp larger ' , ...
                 'than current time. I don''t know how to interpret this']);
            existing = existing{end};
        end

        function printHeader(rt)
            str = sprintf(' Running test %s of group %s ', rt.name, rt.group);
            pad = 8;
            fprintf('\n\n');
            fprintf(repmat('*', 1, numel(str) + 2*pad));
            fprintf('\n');
            fprintf(repmat('*', 1, pad));
            fprintf(str);
            fprintf(repmat('*', 1, pad));
            fprintf('\n');
            fprintf(repmat('*', 1, numel(str) + 2*pad));
            fprintf('\n\n');
        end

        function printFooter(rt, flag)
            str = sprintf(' Test %s %s ', rt.name, rt.flag2str(flag));
            pad = 8;
            fprintf('\n\n');
            fprintf(repmat('*', 1, numel(str) + 2*pad));
            fprintf('\n');
            fprintf(repmat('*', 1, pad));
            fprintf(str);
            fprintf(repmat('*', 1, pad));
            fprintf('\n');
            fprintf(repmat('*', 1, numel(str) + 2*pad));
            fprintf('\n\n');
        end
    end

    methods (Static)
        function fmt = dateTimeFormat()
            fmt = 'yyyy-MM-dd_HH-mm-ss';
        end

        function str = flag2str(flag)
            switch flag
                case -1
                    str = 'inconclusive';
                case 0
                    str = 'failed';
                case 1
                    str = 'passed';
            end
        end
    end
end

function tf = isPackedProblem(x)
    tf = isstruct(x) && isfield(x, 'BaseName') && isfield(x, 'SimulatorSetup') ...
        && isfield(x, 'OutputHandlers');
end

function name = defaultName(setup)
    if ischar(setup) || isstring(setup)
        name = char(setup);
    elseif isa(setup, 'function_handle')
        name = func2str(setup);
    elseif isPackedProblem(setup) && isfield(setup, 'Name')
        name = setup.Name;
    elseif isstruct(setup) && isfield(setup, 'name')
        name = setup.name;
    else
        name = 'regression-test';
    end
end
