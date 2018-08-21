function [F,Phi] = computeFandPhiFromDist(t, vals, varargin)
opt = struct('volume', [], 'allocation', []);
opt = merge_options(opt, varargin{:});
% assure vector shape
[t, vals] = deal(t(:), vals(:));

if vals(end) == 0
    warning('All values of distribution are zero, empty output.')
    [F,Phi] = deal([]);
    return
end

% take mean values of intervals
t(isnan(t)) = 0;
vals(isnan(vals)) = 0;
[mt, mvals] = deal(.5*t(1:end-1)+.5*t(2:end), .5*vals(1:end-1)+.5*vals(2:end));

% integrate to get cum flux
F = cumsum(diff(t).*mvals);
% normalize
if isempty(opt.allocation)
    F = F/F(end);
else
    F = F/opt.allocation;
end

% integrate to get cum volume
Phi = cumsum(diff(t).*mvals.*mt);
% normalize
if isempty(opt.volume)
    Phi = Phi/Phi(end);
else
    Phi = Phi/opt.volume;
end

% Add point (1,1) if not present
if F(end) ~= 1 || Phi(end) ~= 1
    [F, Phi] = deal([F;1], [Phi;1]);
end
end

