function report = comparePackedRegressionResults(rt, problemCurr, problemPrev, details)
%Compare two packed simulation problems using MRST regression rules.

    if nargin < 4 || isempty(details)
        detailParts = {};
    elseif iscell(details)
        detailParts = details;
    else
        detailParts = {details};
    end

    [wsc, stc, repc] = getPackedSimulatorOutput(problemCurr, 'readFromDisk', false);
    [wse, ste, repe] = getPackedSimulatorOutput(problemPrev , 'readFromDisk', false);
    report = struct('passed', true, 'details', nan);

    fun = @(v) norm(v, inf);

    if rt.compareWellSols && ~isempty(wsc) && ~isempty(wsc{1})
        wsRep = compareOutputStructsLocal(wsc, wse, rt.tolerance.wellSols, ...
                                          'fun', fun, 'includeStructs', false);
        if ~wsRep.passed, detailParts{end+1} = 'wellSols'; end %#ok<AGROW>
        report.wellSols = wsRep;
        report.passed = report.passed && wsRep.passed;
    end
    if rt.compareStates
        stRep = compareOutputStructsLocal(stc, ste, rt.tolerance.states, ...
                                         'fun', fun, 'skip', {'wellSol'}, ...
                                         'includeStructs', false);
        if ~stRep.passed, detailParts{end+1} = 'states'; end %#ok<AGROW>
        report.states = stRep;
        report.passed = report.passed && stRep.passed;
    end
    if rt.compareReports
        repRep = compareSimulationReportsLocal(repc, repe, rt.tolerance.reports);
        if ~repRep.passed, detailParts{end+1} = 'reports'; end %#ok<AGROW>
        report.reports = repRep;
        report.passed = report.passed && repRep.passed;
    end

    if isempty(detailParts)
        details = '--';
    else
        details = [strjoin(detailParts, ', '), ...
                   ' differ in magnitude or number of timesteps by more than prescribed tolerance.'];
    end

    report.passed = report.passed*1;
    report.details = details;
end

function report = compareOutputStructsLocal(struct1, struct2, tol, varargin)
    if isa(struct1, 'Resulthandler')
        n1 = struct1.numelData();
        n2 = struct2.numelData();
    else
        n1 = numel(struct1);
        n2 = numel(struct2);
    end

    n = min(n1, n2);
    if n == 0
        report = struct('dsteps', n1 - n2, 'passed', false);
        return
    end

    structTmp = cell(n, 1);
    for i = 1:n
        st1 = struct1{i};
        st2 = struct2{i};
        structTmp{i} = compareStructs(st1, st2, varargin{:});
    end
    structd = struct();
    for name = fieldnames(structTmp{1})'
        v = cellfun(@(st) st.(name{1}), structTmp);
        structd.(name{1}) = v;
    end
    report = struct('dvalues', structd, 'dsteps', n1 - n2);
    report.passed = checkDifferenceLocal(report, tol);
end

function report = compareSimulationReportsLocal(report1, report2, tol)
    n1 = report1.numelData();
    n2 = report2.numelData();
    n = min(n1, n2);
    if n == 0
        report = struct('dsteps', n1 - n2, 'passed', false);
        return
    end
    names = {'nonlinearIterations', 'linearIterations'};
    out = struct();
    for name = names
        out1 = getReportOutput(report1, 'type', name{1});
        out2 = getReportOutput(report2, 'type', name{1});
        dout = out1.total(1:n) - out2.total(1:n);
        out.(name{1}) = dout;
    end
    report = struct('dvalues', out, 'dsteps', n1 - n2);
    report.passed = checkDifferenceLocal(report, tol);
end

function ok = checkDifferenceLocal(report, tol)
    ok = abs(report.dsteps) == 0;
    names = fieldnames(report.dvalues);
    for name = names'
        toltmp = tol;
        if strcmpi(name{1}, 'dvalues')
            toltmp = 0;
        end
        ok = ok & norm(report.dvalues.(name{1}), inf) <= toltmp;
    end
end
