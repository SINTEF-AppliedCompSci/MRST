function r = computePhaseRatesFromRTD(rtd, fluid, varargin)
% Compute time-dependent 2-phase region recovery factor(s) based on RTD
%
% rtd: rtd-structure as computed by computeRTD
% (tau, cums) pairs of tof-values and integrated saturation-fronts at time = 1
% such that at time t we have (t*tau, t*cums).
% We should always have cums(end)=1
%
% If s0 is defaulted, s0=s_conn is assumed, otherwise a crude approximation
% is made as follows:
% Ahead of the front, s0 is kept fixed. At time t, the front is kept fixed,
% but the time t->t* is adjusted such that the the t*-integral of (s(t)-s0)
% equals
%
% RETURNS
% structure r with fields
%   't'  - time vector(s)
%   'rf' - recovery factors cooresponding to t
%   'q'  - total flow rate for each rf(t)
%   'props' - structure containing fluid properties

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

opt = struct('sw0',  [], 'waterCutTol', .999, 'allocationTol', 0.001);
opt = merge_options(opt, varargin{:});

nreg = size(rtd.values, 2);
% check for initial satuation, optional argument has presendence over
% rtd-field, assume WO for now
sw0 = [];
if ~isempty(opt.sw0)
    sw0 = opt.sw0;
elseif isfield(rtd, 'sw0')
    sw0 = rtd.sw0;
end

p_ref = 200*barsa; % fluid should be incompressible here, so pressure immaterial
alloc     = rtd.allocations;
zeroAlloc = alloc/sum(alloc) < opt.allocationTol;      

r = struct('pairIx', rtd.pairIx, 't', rtd.t, 'ratesW', [], 'ratesO', [], ...
           'allocations', rtd.allocations, 'volumes', rtd.volumes, ...
           'volumesW', [], 'volumesO', []);
t    = rtd.t(:,1);
    
% if sw0 is non-uniform, use explicit solver
if numel(sw0) > 1 && rtd.reverse
    assert(size(sw0,2) == size(rtd.values,2), 'Size of initial saturations does not match rtd-struct');
    wcut = computeWatercutFromRTD(rtd, fluid,  'waterCutTol', opt.waterCutTol, ...
        'allocationTol', opt.allocationTol);
else % use Buckley Leverett
    % check format of sw0
    %singleInitSat = true;
    if ~isempty(sw0)
        if numel(sw0)>1
     %       singleInitSat = false;
            if size(sw0, 2) > 1
                disp('Warning: Initial saturation will be averaged')
                sw0 = mean(mean(sw0(:, ~zeroAlloc)));
            end
        end
    end
    [tau, s, fw] = getSaturationFront(fluid, 'sw0', sw0);
    fw = fw(2:end);%.5*(fracFlow(1:end-1)+fracFlow(2:end));
    mtau     = .5*(tau(1:end-1,:)+tau(2:end,:));
   
    mt   = .5*(t(1:end-1) + t(2:end));
    % rtd's are scaled such that integral equals allocation (as fraction)
    % % Accordingly multiply by injector-flux
    q = bsxfun(@times, rtd.values, rtd.injectorFlux(rtd.pairIx(:,1)).');
    
    % use cummulative rates to better preserve allocation when interpolating
    cq = cumsum([zeros(1, size(q,2)); ...
        bsxfun(@times, diff(t), q(2:end,:))]);
    mq = .5*(q(1:end-1,:)+q(2:end,:));
    F    = cell(1, nreg);
    for k = 1:nreg
        F{k} = griddedInterpolant(t, cq(:,k), 'linear', 'none');
    end
    
    wcut = nan(size(rtd.values));
    onlyWater = false(nreg, 1);
    for k = 2:numel(t)
        tcur = t(k);
        %if singleInitSat
            taucur  = tcur*tau;  % current t-scaling of saturation front
            %mtaucur = tcur*mtau;
            mtaucur = tcur*tau(2:end);
            for kr = 1:nreg
                if ~onlyWater(kr) && ~zeroAlloc(kr)
                    cqcur = F{kr}(taucur);
                    %w    = diff(cqcur).*diff(taucur)/alloc(kr);
                    w    = diff(cqcur)/alloc(kr);
                    wcut(k, kr) = sum(fw.*w, 'omitnan');
                    if wcut(k, kr) > opt.waterCutTol
                        onlyWater(kr) = true;
                        wcut(k+1:end, kr) = 1;
                    end
                end
                
            end
        %end
    end
    wcut(1,:) = wcut(2,:);
    wcut = min(wcut, 1, 'includenan');
end
r.ratesW = bsxfun(@times, wcut, alloc.')*fluid.bW(p_ref);
r.ratesW(isnan(r.ratesW)) = 0;
r.ratesO = bsxfun(@times, 1-wcut, alloc.')*fluid.bO(p_ref);
r.ratesO(isnan(r.ratesO)) = 0;
r.allocations(zeroAlloc) = eps;

%% set initial water/oil volumes (at standard conditions)
% sw0 = max(0, sum(.5*(s(1:end-1)+s(2:end)).*diff(tau))-1)/tau(end);
r.volumesW = rtd.volumesW*fluid.bW(p_ref);
r.volumesO = (r.volumes-r.volumesW)*fluid.bO(p_ref);


end
