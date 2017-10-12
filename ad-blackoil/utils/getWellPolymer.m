function [wPoly, wciPoly, iInxW] = getWellPolymer(W)
    if isempty(W)
        wPoly = [];
        wciPoly = [];
        iInxW = [];
        return
    end
    inj = vertcat(W.sign) == 1;
    injix   = find(inj);
    polInj = cellfun(@(x)~isempty(x), {W(injix).c});
    wPoly = zeros(nnz(injix), 1);
    wPoly(polInj) = vertcat(W(injix(polInj)).c);
    wciPoly = rldecode(wPoly, cellfun(@numel, {W(injix).cells}));

    % Injection cells
    nPerf = cellfun(@numel, {W.cells})';
    nw    = numel(W);
    perf2well = rldecode((1:nw)', nPerf);
    compi = vertcat(W.compi);
    iInx  = rldecode(inj, nPerf);
    iInx  = find(iInx);
    iInxW = iInx(compi(perf2well(iInx),1)==1);
end

