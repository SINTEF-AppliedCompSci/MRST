function [F,Phi] = computeFandPhiFromDist(RTD, varargin)
opt = struct('sum',false);
opt = merge_options(opt, varargin{:});

if strcmpi(RTD.creator,'estimateRTD')
    warning(['F-Phi diagrams computed from estimated RTD ', ...
        'are highly inaccurate. Proceed with caution\n']);
    normalizeAlloc = true;
    normalizeVol = false;
else
    normalizeAlloc = false;
    normalizeVol = true;
end

if opt.sum
    t     = RTD.t(:,1);
    vals  = sum(RTD.values,2);
    vol   = sum(RTD.volumes);
    alloc = sum(RTD.allocations);
else
    t     =  RTD.t;
    vals  = RTD.values;
    vol   = RTD.volumes(:).';
    alloc = RTD.allocations(:).';
end

% take mean values of intervals
t(isnan(t)) = 0;
vals(isnan(vals)) = 0;
[mt, mvals] = deal(.5*t(1:end-1,:)+.5*t(2:end,:),vals(2:end,:));

% integrate to get cum flux
F = cumsum(diff(t,1).*mvals, 1);
if normalizeAlloc
    F = bsxfun(@rdivide, F, F(end,:));
else
    F = bsxfun(@rdivide, F, alloc);
end

% integrate to get cum volume
Phi = cumsum(diff(t,1).*mvals.*mt,1);
if normalizeVol
    Phi = bsxfun(@rdivide, Phi, vol);
else
    Phi = bsxfun(@rdivide, Phi, Phi(end,:));
end

% Add points (0,0) and (1,1) to be sure they exist
one  = ones (1,numel(vol));
zero = 0*one;
[F, Phi] = deal([zero; F; one], [zero; Phi; one]);
end