classdef RegressionTest
    
    properties
        group = 'general'
        deleteExisting = true
        name
        test
        compareFn = []
        tolerance = 0
        testCaseOpt = {}
        problemInput = {}
        verbose = true
    end
    
    methods
        %-----------------------------------------------------------------%
        function rt = RegressionTest(test, varargin)
            [rt, testCaseOpt] = merge_options(rt, varargin{:});

            if ischar(test)
                name = test;
            elseif isa(test, 'TestCase')
                name = test.name;
            else
                error(['Input argument `test` must be either a '     , ...
                       'string corresponding to a test case setup '  , ...
                       'function, or an object of the TestCase class']);
            end
            if isempty(rt.name), rt.name = name; end
            rt.test = test;
            rt.testCaseOpt = testCaseOpt;
            
            if isempty(rt.compareFn)
                rt.compareFn = @(test, pc, pe) compareResults(test, pc, pe, opt);
            end
        end
    
        %-----------------------------------------------------------------%
        function report = runRegressionTest(rt)
            time = rt.getTimeString();
            testCase = rt.test;
            if ischar(testCase)
                testCase = TestCase(rt.test, rt.testCaseOpt{:});
            end
            if rt.verbose, rt.printHeader(); end
            % Output directory for test case
            directory = rt.getDataDirectory();
            % Get test case hash
            hash = rt.getRegressionTestHash(testCase);
            % Find existing resutls for test case
            nameExt = rt.findExisting(hash, time);
            if ~isempty(nameExt)
            % Existing resuls found - get result handler
            probExt = testCase.getPackedSimulationProblem(                 ...
                              'Directory', directory                       , ...
                              'Name'     , nameExt, rt.problemInput{:});
            else
                % No exising results found. Issue warning saying that this run
                % will be inconclusive
                warning(['No existing results for test case %s. Results '    , ...
                         'of this test will be inconclusive.'], testCase.name);
            end
            % Get name for current result folder
            nameCurr = [time, '_', hash];
            % Set up packed problem for current results
            probCurr  = testCase.getPackedSimulationProblem(             ...
                              'Directory', directory                     , ...
                              'Name'     , nameCurr, rt.problemInput{:});
            % Simulate
            try
                simulatePackedProblem(probCurr);
                details = [];
            catch ex
                details = ex.message;
            end
            % Compare to existing results if they exist
            if ~isempty(nameExt)
                report = rt.compareResults(probCurr, probExt, details);
                if report.passed == 1 && rt.deleteExisting
                    % Delete existing results
                    rmdir(fullfile(directory, nameExt), 's');
                end
            else
                if ~isempty(details), details = [details, '. ']; end
                details = [details, 'No existing results to compare againts'];
                report  = struct('passed', -1, 'details', details);
            end
            if rt.verbose, rt.printFooter(report.passed); end
            % Save regression test report
            save(fullfile(directory, ['regreport_', nameCurr, '.mat']), 'report');
        end
        
        %-----------------------------------------------------------------%
        function reports = getTestReports(rt)
            directory = rt.getDataDirectory();
            l = ls(directory);
            l = regexp(char(l), '\s+', 'split');
            reports = cell(numel(l), 1);
            for i = 1:numel(l)
                if ~contains(l{i}, 'regreport_'), continue; end
                data = load(fullfile(directory, l{i}));
                reports{i} = data.report;
            end
            reports = reports(~cellfun(@isempty, reports));
        end
        
        %-----------------------------------------------------------------%
        function hash = getRegressionTestHash(rt, testCase)
            hash   = testCase.getTestCaseHash();
            if isempty(rt.problemInput), return; end
            piHash = cellfun(@(ip) obj2hash(ip), rt.problemInput, 'UniformOutput', false);
            hash   = str2hash(strjoin([hash, piHash], '_'));
        end
        
        %-----------------------------------------------------------------%
        function directory = getDataDirectory(rt)
            directory = fullfile(mrstOutputDirectory(), 'reg-tests', rt.group(), rt.name);
        end
        
        %-------------------------------------------------------------------------%
        function time = getTimeString(rt)
            time = char(datetime('now', 'Format', rt.dateTimeFormat()));
        end
        
        %-----------------------------------------------------------------%
        function existing = findExisting(rt, hash, timeCurr)
            
            directory = rt.getDataDirectory();
            
            existing = [];
            if ~exist(directory, 'dir'); return; end
            dirNames = ls(directory);
            if isempty(dirNames), existing = []; return; end

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
            timeExt  = cellfun(@(ext) datetime(ext(1:numel(fmt)), ...
                                                'Format', fmt), existing);
            if numel(timeExt) > 1
                warning(['More than one set of existing results found. '    , ...
                         'Using most recent results for regression testing.']);

                [timeExt, ix] = sort(timeExt);
                existing      = existing(ix);

                prompt = sprintf('Do you want to delete older results? y/n [y]: ');
                str    = input(prompt,'s');
                if strcmpi(str, 'n')
                    fprintf('Ok, will not remove results.\n');
                else
                    for i = 1:numel(existing)-1
                        fn = fullfile(directory, existing{i});
                        rmdir(fn, 's');
                    end
                end
            end
            timeExt = timeExt(end);
            assert(timeExt(end) < timeCurr, ...
                   ['Previous results are stored with timestamp larger '   , ...
                    'than current time. I don''t know how o interpret this']);
            existing = existing{end};
        end
        
        %-----------------------------------------------------------------%
        function report = compareResults(rt, problemCurr, problemExt, details)
            % Parse deails
            if ~isempty(details), details = [details, '. ']; end
            % Get result handlers
            [wsc, stc, repc] = getPackedSimulatorOutput(problemCurr, 'readFromDisk', false);
            [wse, ste, repe] = getPackedSimulatorOutput(problemExt , 'readFromDisk', false);
            % Use inf norm for comparing all fileds
            fun = @(v) norm(v, inf);
            % Compare states
            stRep = rt.compareStructsLocal(stc, ste, rt.tolerance,      ...
                                         'fun'           , fun        , ...
                                         'omit'          , {'wellSol'}, ...
                                         'includeStructs', false      );
            if ~stRep.passed, details = [details, 'states']; end
            % Compare well solutions
            wsRep = rt.compareStructsLocal(wsc, wse, rt.tolerance,  ...
                                           'fun'           , fun  , ...
                                           'includeStructs', false);
            if ~wsRep.passed, details = [details, {'wellSols'}]; end
            % Compare simulation reports
            repRep = rt.compareReports(repc, repe, 0);
            if ~repRep.passed, details = [details, {'reports'}]; end
            if isempty(details)
                details = '--';
            else
                details = strjoin(details, ', ');
                details = [details, ' differ in magnitude or number of '         , ...
                                    'timesteps by more than prescribed tolerance.'];
            end
            % Make report
            passed = (stRep.passed && wsRep.passed && repRep.passed)*1;
            report = struct('passed'  , passed , ...
                            'details' , details, ...
                            'states'  , stRep  , ...
                            'wellSols', wsRep  , ...
                            'reports' , repRep);
        end
        
        %-----------------------------------------------------------------%
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
        
        %-----------------------------------------------------------------%
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
        %-----------------------------------------------------------------%
        function fmt = dateTimeFormat()
            fmt = 'yyyy-MM-dd_HH-mm-ss';
        end
        
        %-----------------------------------------------------------------%
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
        
        %-----------------------------------------------------------------%
        function report = compareStructsLocal(struct1, struct2, tol, varargin)
            n1 = struct1.numelData();
            n2 = struct2.numelData();

            n  = min(n1, n2);
            if n == 0, report = struct('dsteps', n1 - n2, 'passed', false); return; end
            structd = [];
            for i = 1:n
                st1 = struct1{i};
                st2 = struct2{i};
                structTmp = compareStructs(st1, st2, varargin{:});
                structd   = addStates(structd, structTmp);
            end
            report = struct('dvalues', structd, 'dsteps', n1 - n2);
            report.passed = checkDifference(report, tol);
        end

        %-----------------------------------------------------------------%
        function report = compareReports(report1, report2, tol)
            n1 = report1.numelData();
            n2 = report2.numelData();
            if min(n1, n2) == 0, report = struct('dsteps', n1 - n2, 'passed', false); return; end
            names = {'nonlinearIterations', 'linearIterations'};
            out = struct();
            for name = names
                out1 = getReportOutput(report1, 'type', name{1});
                out2 = getReportOutput(report2, 'type', name{1});
                n = min(numel(out1.total), numel(out2.total));
                dout = norm(out1.total(1:n) - out2.total(1:n), inf);
                out.(name{1}) = dout;
            end
            report = struct('dvalues', out, 'dsteps', n1 - n2);
            report.passed = checkDifference(report, tol);
        end
        
    end
        
end

% Helpers
 %-----------------------------------------------------------------%
function ok = checkDifference(report, tol)
    ok = abs(report.dsteps) == 0;
    names = fieldnames(report.dvalues);
    for name = names'
        ok = ok & report.dvalues.(name{1}) <= tol;
    end
end

%-----------------------------------------------------------------%
function st = addStates(st1, st2)
    st = st2;
    if isempty(st1), return; end
    fnames = fieldnames(st1);
    for fname = fnames'
        st.(fname{1}) = st1.(fname{1}) + st2.(fname{1});
    end
end