function alpha = calcAlpha(W, bq, b, rs, perf2well)
% finally explicit computations of alpha:

if nargin < 5 % single well
    perf2well = ones(numel(W.cells), 1);
end

qs   = [double(bq{1}), double(bq{2}), double(bq{3}) + double(rs).*double(bq{2})];
rs   = double(rs);
b    = [double(b{1}), double(b{2}), double(b{3})];

alpha = cell(numel(W),1);
for wnr = 1:numel(W);
    w = W(wnr);

    pinx = (perf2well == wnr);
    qsw   = qs(pinx,:);
    rsw   = rs(pinx);
    bw    = b(pinx,:);

    nperf = numel(w.cells);
    ii = [w.topo(:,2); w.topo(2:end, 1)];
    jj = [(1:nperf)'; (2:nperf)'];
    vv = [ones(nperf, 1); -ones(nperf-1, 1)];
    C = sparse(ii, jj, vv, nperf, nperf);

    wbqs = C\qsw;
    % take flow in each segment as equal to upstream:
    % reindex to avoid 0-index (producers)
    segInflows = zeros(nperf,3);
    for ph = 1:3
        dwnFlow   = wbqs(:,ph)>0;
        segInflows(w.topo(dwnFlow,2),ph) =  wbqs(dwnFlow,ph);

        upFlow    = [false; wbqs(2:end,ph)<0];
        segInflows(w.topo(upFlow,1),ph)  = -wbqs(upFlow,ph);

        isPrd  = qsw(:,ph)<0;
        segInflows(isPrd, ph) = segInflows(isPrd, ph) - qsw(isPrd,ph);
    end
    % at connection conditions:
    segq   = [segInflows(:,1)./bw(:,1), ...
              segInflows(:,2)./bw(:,2), ...
             (segInflows(:,3) - rsw.*segInflows(:,2))./bw(:,3)];
    % alpha{wnr}  = abs(segq)./(sum(abs(segq),2)*[1 1 1]);

    alpha{wnr} = repmat([1 0 0], nperf, 1);
    significant_flux = sum(abs(segq),2) > 1e-6;  % is this small enough?
    alpha{wnr}(significant_flux, :)  = abs(segq(significant_flux, :))./(sum(abs(segq(significant_flux, ...
                                                      :)),2)*[1 1 1]);

end
end