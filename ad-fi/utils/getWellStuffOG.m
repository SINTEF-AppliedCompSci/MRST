function [Tw, dzw, Rw, wc, perf2well, pInx, iInxG, iInxO] = getWellStuffOG(W)
if isempty(W)
    [Tw, dzw, Rw, wc, perf2well, pInx, iInxG] = deal([]);
    return
end

nPerf = cellfun(@numel, {W.cells})';
nPerfTot = sum(nPerf);
nw    = numel(W);
perf2well = rldecode((1:nw)', nPerf);

Rw = sparse((1:nPerfTot)', rldecode((1:nw)', nPerf), 1, nPerfTot, nw);

%------------------------------------------

Tw    = vertcat(W.WI);
dzw   = vertcat(W.dZ);
wc    = vertcat(W.cells);
inj   = vertcat(W.sign)==1;
compi = vertcat(W.compi);
iInx  = rldecode(inj, nPerf);
pInx  = ~iInx;
iInx   = find(iInx);
pInx   = find(pInx);
iInxG  = iInx(compi(perf2well(iInx),3)==1);
iInxO  = iInx(compi(perf2well(iInx),2)==1);
%------------------------------

end
