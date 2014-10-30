function [qWs, qOs, qGs, bhp] = wellSolToVector(wellsols)
% Helper function which makes cell arrays of well solutions easier to plot
    nt = numel(wellsols);
    nw = numel(wellsols{1});
    ws = vertcat(wellsols{:});

    fix = @(v) reshape(v, nw, nt)';
    sgn = fix([ws.sign]);
    bhp = fix([ws.bhp]);
    qWs = sgn.*fix([ws.qWs]);
    qGs = sgn.*fix([ws.qGs]);
    qOs = sgn.*fix([ws.qOs]);

end
