function [withinLims, W] = checkLims(W, pBH, qWs, qOs, qGs)
nwells = numel(W);
withinLims = true(nwells,1);

pBH = double(pBH);
qs = [double(qWs), double(qOs), double(qGs)];


for wnr = 1:numel(W)
    w = W(wnr);
    pBHw = pBH(wnr);
    qsw  = qs(wnr,:);

    if ~isnumeric(w.lims)
        if w.sign > 0   % injector
            modes   = {'bhp', 'rate'};
            flags = [pBHw > w.lims.bhp, sum(qsw) > w.lims.rate];
        else            % producer
            modes   = {'bhp', 'orat', 'lrat', 'grat', 'wrat'};
            flags = [pBHw < w.lims.bhp,           ...
                qsw(2) < w.lims.orat,        ...
                qsw(1)+qsw(2) < w.lims.lrat,   ...
                qsw(3) < w.lims.grat,        ...
                qsw(1) < w.lims.wrat];
        end
    else
        modes = {};
        flags = false(numel(W), 1);
        assert(isinf(w.lims))
    end
    %limits we need to check (all others than w.type):
    chkInx = ~strcmp(w.type, modes);
    vltInx = find(flags(chkInx), 1);
    if ~isempty(vltInx)
        withinLims = false;
        modes  = modes(chkInx);
        switchMode = modes{vltInx};
        fprintf('Well %s: Control mode changed from %s to %s.\n', w.name, w.type, switchMode);
        w.type = switchMode;
        w.val  = w.lims.(switchMode);
        W(wnr) = w;
    end
end
end
