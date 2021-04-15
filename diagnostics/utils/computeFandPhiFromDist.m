function [F,Phi] = computeFandPhiFromDist(RTD, varargin)
%Undocumented Utility Function

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

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
    t       = RTD.t(:,1);
    vals    = sum(RTD.values,2);
    vol     = sum(RTD.volumes);
    alloc   = sum(RTD.allocations);
    injflux = sum(RTD.injectorFlux(:));
else
    t       = RTD.t;
    vals    = RTD.values;
    vol     = RTD.volumes(:).';
    alloc   = RTD.allocations(:).';
    injflux = RTD.injectorFlux(RTD.pairIx(:,1)).';
end

% take mean values of intervals
t(isnan(t)) = 0;
vals(isnan(vals)) = 0;
[mt, mvals] = deal(.5*t(1:end-1,:)+.5*t(2:end,:),vals(2:end,:));

% Integral of RTD should approximate fractional recovery (produced mass/injected mass) 
% Cummulative allocation flux is integral of RTD times injector flux
F = cumsum(diff(t,1).*mvals, 1);
F = bsxfun(@times, F, injflux);
if normalizeAlloc
    F = bsxfun(@rdivide, F, F(end,:));
else
    F = bsxfun(@rdivide, F, alloc);
end

% interaction volume = interaction allocation x mean TOF
% integrate and multiply by allocation to get cummulative volume
Phi = cumsum(diff(t,1).*mvals.*mt,1);
Phi = bsxfun(@times, Phi, alloc);
if normalizeVol
    Phi = bsxfun(@rdivide, Phi, Phi(end,:));
else
    Phi = bsxfun(@rdivide, Phi, vol);
end

% Add points (0,0) and (1,1) to be sure they exist
one  = ones (1,numel(vol));
zero = 0*one;
[F, Phi] = deal([zero; F; one], [zero; Phi; one]);
end
