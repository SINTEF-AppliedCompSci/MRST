function [wcuts, allocations] = computeWatercutFromRTD(rtd, fluid, varargin)
% Compute region watercut for produced fluid by solving 1D transport
% problem along oppsite direction of producer backward TOF.  
%
% SYNOPSIS:
%  [wcuts, allocations] = computeWatercutFromRTD(rtd, fluid, [ix])
%
% DESCRIPTION:
% For a (set of) interaction region(s), this function solves the 1D transport
% problem and outputs watercuts corresponding to time-steps corresponding to 
% rtd.t . The time-steps are assumed to be produced by computeRTD such that the
% step sequence is of the form 
% [dtau0*ones(1,n), m*dtau0*ones(1,n), m^2*dtau0*ones(1,n), ..., m^(N-1)*dtau0*ones(1,n_end)]
% with preferably mod(n,m)=0. 
% For the N periods above the satuation is solved on a gradually coarsened tau-grid, 
% such that for period k, the smallest tau-interval equals m^(k-1)*dtau0
% and hence a stable time-step is m^(k-1)*dtau0/max(dfw).
%
% REQUIRED PARAMETERS:
% rtd   - RTD-structure as copmuted by computeRTD
% fluid - Preferably incompressible water-oil fluid structure
%
% OPTIONAL PARAMETERS
% ix    - index to subset of rtd's. Default all.
%
% 
% RETURNS:
% wcut -       matrix where column k is the estimated wcuts for rtd number ix(k) 
%              at times corresponding to rtd.t
% allocations- corresponding allocations

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

opt = struct('pairIx',                    [], ...
             'allocationTol',     1e-3, ...
             'waterCutTol', .995, ...
             'referencePressure',  200*barsa);        
opt = merge_options(opt, varargin{:});

assert(isfield(rtd, 'sw0')&&rtd.reverse, ...
       'Function computeWatercutFromRTD requires rtd including satuaration along backward TOFS (rtd.swr)');

if isempty(opt.pairIx)
    opt.pairIx = 1:size(rtd.values,2);
end
ix = opt.pairIx;

% Get list of period time-steps and number of steps in each period
[dts, nsteps] = getSimSettings(rtd);

% get fracflow function
pr = opt.referencePressure;  % typically redundant here (constant viscosities)
[muw, muo] = deal(fluid.muW(pr), fluid.muO(pr)); 
if isfield(fluid, 'krO')
    krO = fluid.krO;
else
    krO = fluid.krOW;
end
[lw, lo]   = deal(@(s)fluid.krW(s)/muw, @(s)krO(1-s)/muo);
fw = @(s)lw(s)./(lw(s)+lo(s));

% find max derivative of fw
sw_tmp = initVariablesAD_diagonal((0:.001:1)');
fw_tmp = fw(sw_tmp);
dfwMax = max(fw_tmp.jac{1}.diagonal);
clear sw_tmp fw_tmp

% resample and and setup interpolant for speedup
sw_tmp = (0:.001:1)';
fw = griddedInterpolant(sw_tmp, fw(sw_tmp));

% preallocate
wcuts       = nan( numel(rtd.t(:,1)), numel(ix) );
allocations = rtd.allocations(ix);
% approximate number of saturation steps (trim afterwards):
[ts, vals] = deal( nan( sum( ceil(nsteps*dfwMax)+1), 1) );

% solve a 1D transport problem along reverse TOF with source given by RTD 
t        = rtd.t(:, 1);
totAlloc = sum(allocations);
for reg = 1:numel(ix)
    if allocations(reg) >= opt.allocationTol*totAlloc
        sw = rtd.sw0(1:end-1,ix(reg));
        % scale rtd such that integral to infinity equals 1
        q   = rtd.values(:,ix(reg))*( (rtd.injectorFlux(rtd.pairIx(ix(reg),1))/rtd.allocations(ix(reg))) );
        % don't solve for possible zero fluxes at high tofs 
        [q, sw, tcur, dtscur, nstepscur] = reduceToRelevant(q, sw, t, dts, nsteps);
        
        src = q(2:end).*diff(tcur);
        % flux (rate) along tof equals integral of src. Adjust for "missing" source due to
        % finite time
        adjust   = 1-sum(src);
        src(end) = src(end) + adjust;
        flux = max(0, cumsum( src , 'reverse'));
       
        [ts(1), vals(1)] = deal(0, fw(sw(1)));
        cnt = 0;
        for k_period = 1:numel(dtscur)
            dt = dtscur(k_period);
            % resample (upscale) such that min(diff(t))>=dt
            [tcur, sw, flux, src] = resample(tcur, sw, flux, src, dt);
            % influx plays the role as "porosity"
            poro = [flux(2:end); src(end)];
            curDT = dt/dfwMax; % stable time step
            curT  = 0;
            curTTot = dt*nstepscur(k_period);
            while curT < curTTot
                step  = min(curDT, curTTot-curT);
                curT  = curT + step;
                fluxW = [fw(sw).*flux; 0];
                inc   = (src + diff(fluxW))./diff(tcur);
                sw    = sw + step*(inc./poro); 
                %plot(tcur(1:end-1)/year, sw)
                %drawnow
                cnt = cnt+1;
                wc = fw(sw(1));
                [ts(cnt+1), vals(cnt+1)] = deal(ts(cnt)+step, wc);
                if wc > opt.waterCutTol
                    break
                end
            end
        end
        [ts, vals] = deal(ts(1:(cnt+1)), vals(1:(cnt+1)));
        % interpolate onto rtd t-values
        wcuts(:, reg) = interp1(ts, vals, rtd.t(:,1), 'linear', 1);
    end
end
end

% -------------------------------------------------------------------------
function [dts, nsteps] = getSimSettings(rtd)

dt = diff(rtd.t(:,1));
[dts, ~, ix] = uniquetol(dt);
[~, nsteps]  = rlencode(ix);
% check
if numel(dts) > 1
    mlt = dts(2:end)./dts(1:end-1);
    assert(all(mlt>1) && numel(uniquetol(mlt))==1, ...
           'RTD t-values of unexpected type');
end
end

% -------------------------------------------------------------------------
function [t, sw, flux, src] = resample(t, sw, flux, src, dtMin)
    % resample (upscale) such that min(diff(t))>=dtMin
    dt = diff(t);
    if dtMin/min(dt) <= 1+sqrt(eps)
        % do nothing
    else
        ix = nan(numel(dt),1);
        bin = 1;
        sumt = 0;
        for k = 1:numel(dt)
            ix(k) = bin; 
            sumt = sumt + dt(k);
            if sumt/dtMin >= 1-sqrt(eps)
                bin  = bin+1;
                sumt = 0;
            end
        end
        % dt-weghted saturation average 
        %sw = accumarray(ix, sw.*dt)./accumarray(ix, dt);
        % sum sources
        src = accumarray(ix, src);
        %poro = [flux(2:end); src(end)];
        %sw = accumarray(ix, sw.*poro.*dt)./accumarray(ix, poro.*dt);
        [~, n] = rlencode(ix);
        lastIx = cumsum(n);
        sw = [sw(1); sw(lastIx(1:end-1)+1)];
        t = [t(1); t(lastIx+1)];
        % flux at new t-values 
        flux = [flux(1); flux(lastIx(1:end-1)+1)];
    end
end
% -------------------------------------------------------------------------
function [q, sw, t, dts, nsteps] = reduceToRelevant(q, sw, t, dts, nsteps)
ix = find( cumsum(q(2:end).*diff(t))<= 1-1e-7, 1, 'last' );
if ix < numel(q)-1
    %[tcur, sw, flux, src, dtscur, nstepscur] = reduceToRelevant(ix, t, sw, flux, src, dts, nsteps)
    [q, sw, t] = deal(q(1:(ix+1)), sw(1:ix), t(1:(ix+1)));
    outerStepIx = find(cumsum(nsteps) > ix, 1, 'first');
    dts    = dts(1:outerStepIx);
    nsteps = nsteps(1:outerStepIx);
     nsteps(end) = ix - sum(nsteps(1:(outerStepIx-1)));
end
end
