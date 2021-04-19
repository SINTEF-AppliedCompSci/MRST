function dist = computeRTD(state, G, pv, D, WP, W, varargin)
% Compute residence time distributions based on computed tof-values
% SYNOPSIS:
%   dist = computeRTD(state, G, pv, D, WP, W, 'pn1', pv1, ...)
% 
% DESCRIPTION:
%  This function computes well-pair RTDs by simulating an imaginary tracer 
%  which ditributes equally among all phases (follows the total flux-field), 
%  and follows the flux field given by state. A tracer pulse is injected at 
%  injector i at time zero and produced at producer p. The RTD is scaled 
%  such that
%     RTD_{ip}(t) = (rate producer p / total tracer mass) * c_p(t), 
%  where c_p(t) is the tracer concentration in producer p. 
% 
%  With the above scaling the RTD has units [s]^-1
%   * the integral of the RTD approximates fractional recovery 
%    (produced mass/injected mass) 
%   * the mean of the RTD time i-to-p allocation approximates the i-to-p 
%     allocation volume
%
% RETURNS:
% dist - structure with fields
%    pairIx            nreg x 2 each row index to injector/producer 
%    t                 nbin x 1 discrete times
%    values            nbin x 1 discrete RTD-values
%    volumes           nreg x 1 interaction volume for each well pair
%    allocations       nreg x 1 interaction allocation for each well pair
%    injectorFlux      ninj injector total rates
%    producerFlux      nprod producer total rates
%
% SEE ALSO
%  estimateRTD
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


opt = struct('injectorIx', [], ...
             'producerIx', [], ...
             'nsteps',     50, ...
             'nbase',       5, ...
             'computeSaturations', true, ...
             'wbHandle',   [], ...
             'showWaitbar',  true, ...
             'reverse',    false);

opt = merge_options(opt, varargin{:});

[pix, iix] = deal(opt.producerIx, opt.injectorIx);         
if isempty(pix), pix = (1:numel(D.prod))'; end
if isempty(iix), iix = (1:numel(D.inj))'; end       

if opt.nsteps < opt.nbase
    opt.nsteps = opt.nbase;
    fprintf('nsteps-option set equal to nbase-option for convenience ...');
end

% create output struct
nreg = numel(pix).*numel(iix);
dist = struct('pairIx',         nan(nreg, 2),       ...
              't',              [],                 ...
              'volumes',        nan(nreg, 1),       ...
              'allocations',    nan(nreg,1),        ...
              'values',         [],                 ...
              'injectorFlux',   nan(numel(iix),1),  ...
              'producerFlux',   nan(numel(pix),1)); 
          
% get well rates from allocations
for ik = 1:numel(iix)
    dist.injectorFlux(ik) = sum( sum(WP.inj(iix(ik)).alloc) );
end

for pk = 1:numel(pix)
    dist.producerFlux(pk) = sum( sum(WP.prod(pix(pk)).alloc) );
end

% collect relevant data from WP
ix = 0;          
for ik = 1:numel(iix)
    for pk = 1:numel(pix)
        ix = ix +1;
        dist.pairIx(ix,:) = [iix(ik), pix(pk)];
        ixWP = WP.pairIx(:,1)==iix(ik) & WP.pairIx(:,2) == pix(pk);
        dist.volumes(ix)      = WP.vols(ixWP);
        dist.allocations(ix)  = sum(WP.inj(iix(ik)).alloc(:,pix(pk)));
    end
end

% index into well-struct
[inj, prod] = deal(D.inj(iix), D.prod(pix));

% extract relevat subset
sub = sum(D.itracer(:, iix), 2) > 1e-6;

% divide into two or more simulation periods
% 1. 0  <= t <   t1 = min(vol/alloc)  , i.e., minimal interaction mean RT
% 2. t1 <= t < (nbase+1)*t1
%     |
% n. (nbase+1)^(n-2)*t1 <= t <  ~maxtof

% only check allocations > total*10^-3
checkIx = dist.allocations > 1e-3*sum(dist.allocations);

if ~checkIx
    return
end
% end of first period
t1 = min(dist.volumes(checkIx)./dist.allocations(checkIx));
% total pvi
pvi = sum(dist.volumes)/sum(dist.allocations);

if isfield(D, 'itof')
    itof =  D.itof(sub, iix);
else
    itof = D.tof(sub,1);
end
itof(~isfinite(itof)) = 0;
maxTOF = max(max(itof));
% end time
t_end = min(maxTOF, 100*pvi);

% find number of simulation periods and period dt and nsteps
nperiods    = ceil( log( t_end*(opt.nbase-1)/t1+1 )/(log(opt.nbase)) );
dts         = (t1/opt.nsteps)*opt.nbase.^((0:(nperiods-1))');
nsteps      = opt.nsteps*ones(nperiods,1);
nsteps(end) = ceil( (t_end-sum(dts(1:end-1).*nsteps(1:end-1)))/dts(end) );
nStepsTot   = sum(nsteps);

% allocate 
[dist.t, dist.values] = deal(nan(nStepsTot+1, nreg));

% setup system
[sysmat, qp_well, tr0, pv, cix] = setupSystemComponents(state, G, pv, W, sub, inj, prod, opt.reverse);

tr = tr0;
[ni, np] = deal(numel(iix), numel(pix));
if opt.reverse
    [ni, np] = deal(np, ni);
end
vals = cell(1, ni);
for k = 1:numel(vals)
    vals{k} = zeros(nStepsTot+1, np);
end

if opt.computeSaturations
    % get region total and water volumes
    if ~opt.reverse
        pv_reg  = bsxfun(@times, D.ptracer(cix,pix), pv);
        pvw_reg = bsxfun(@times, pv_reg, state.s(cix,1));
    else
        pv_reg  = bsxfun(@times, D.itracer(cix,pix), pv);
        pvw_reg = bsxfun(@times, pv_reg, state.s(cix,1));
    end
    % initial saturation
    sats = vals; 
    for tn = 1:numel(sats)
        curvals  = (tr0(:, tn)'*pvw_reg)./(tr0(:, tn)'*pv_reg);
        curvals(~isfinite(curvals)) = 0;
        sats{tn}(1,:) = curvals;
    end
end

if opt.showWaitbar
    if isempty (opt.wbHandle)
        h = waitbar(0, 'Computing distribution(s) ...');
    else
        h = waitbar(0, opt.wbHandle, 'Computing distribution(s) ...');
    end
end
% Main loop ---------------------------------------------------------------
cnt = 0;
for ti = 1:numel(dts)
    A = sysmat(dts(ti));
    for k = 1:nsteps(ti)
        cnt = cnt+1;
        tr = A\tr;
        for tn = 1:numel(vals)
            curvals  = -(tr(:, tn)'*qp_well);
            curvals(~isfinite(curvals)) = 0;
            vals{tn}(cnt+1,:) = curvals;
            if opt.computeSaturations
                curvals  = (tr(:, tn)'*pvw_reg)./(tr(:, tn)'*pv_reg);
                curvals(~isfinite(curvals)) = 0;
                sats{tn}(cnt+1,:) = curvals;
            end
            if opt.showWaitbar
                waitbar(cnt/((nperiods+1)*opt.nsteps), h);
            end
        end
    end
end
if isempty(opt.wbHandle)&&opt.showWaitbar, close(h); end

% rearrange distributions and possibly scale such that distribution integral
% equals allocation/injectorflux
dt = rldecode(dts,nsteps);
t = cumsum([0; vertcat(dt)]);
dist.t = repmat(t, [1, nreg]);
ix = 0;          
for ik = 1:numel(iix)
    for pk = 1:numel(pix)
        ix = ix +1;
        if ~opt.reverse
            dist.values(:, ix)    = vals{ik}(:, pk);
        else
            fac = -dist.producerFlux(dist.pairIx(ix,2))/dist.injectorFlux(dist.pairIx(ix,1));
            dist.values(:, ix)    = vals{pk}(:, ik)*fac;
        end
    end
end

% rearrage saturation values and store interaction region water volumes
if opt.computeSaturations
    dist.sw0      = nan(nStepsTot+1, nreg);   
    dist.volumesW = nan(size(dist.volumes));
    ix = 0;
    for ik = 1:numel(iix)
        for pk = 1:numel(pix)
            ix = ix +1;
            V_ip  = D.itracer(cix,ik).*D.ptracer(cix,pk).*pv;
            Vw_ip = V_ip.*state.s(cix, 1);
            if ~opt.reverse
                dist.sw0(:, ix)    = sats{ik}(:, pk);
            else
                dist.sw0(:, ix)    = sats{pk}(:, ik);
            end
            dist.volumesW(ix) = sum(Vw_ip);
        end
    end
end
dist.creator = mfilename;
dist.reverse = opt.reverse;
end

% ------------------------------------------------------------------------
function [sysmat, qp_well, tr0, pv, cix] = setupSystemComponents(state, G, pv, W, sub, inj, prod, reverse)
if reverse
    [inj, prod] = deal(prod, inj);
end
nc = G.cells.num;
N  = G.faces.neighbors;
if size(N, 1) ~= size(state.flux, 1) % flux on interior faces
    N = N(~any(N==0, 2),:);
end
% get index to relevant faces
sub = [false; sub];
fix = sub(N(:,1)+1).*sub(N(:,2)+1) > 0;
% redfine N to cix
N   = N(fix,:);
% find corresponding cells (possibly different than sub)
cix = false(nc,1);
cix(N(:)) = true;

% remap to subset
remap = zeros(nc,1);
remap(cix) = (1:nnz(cix))';
N = remap(N);
% swap to get upstream in neighbor-list
if ~reverse
    v =  sum(state.flux(fix,:),2);
else
    v = -sum(state.flux(fix,:),2);
end
neg = v < 0;
N(neg,:) = N(neg, [2 1]);
v(neg)   = -v(neg);

% setup production sources
if ~reverse
    q =  sum(vertcat(state.wellSol.flux), 2);
else
    q = -sum(vertcat(state.wellSol.flux), 2);
end
wcells  = vertcat((W.cells));
qp = sparse(wcells, 1, q.*(q<0), G.cells.num, 1);
qp = qp(cix);

% system matrix
ncr = nnz(cix);
A = sparse(N(:,2), N(:,1), v, ncr, ncr);
d = sum(A, 1)' - qp;
A = A - spdiags(d, 0, ncr, ncr);
A = spdiags(1./pv(cix), 0, ncr, ncr)*A;
% setup system matrix as function of dt
if ~verLessThan('matlab', '9.4')
    sysmat = @(dt)decomposition(speye(ncr)-dt*A);
else
    sysmat = @(dt)(speye(ncr)-dt*A);
end
% setup relevant well-sources columnwise 
[qi_well, qp_well] = deal(sparse(ncr, numel(inj)), sparse(ncr, numel(prod)));
for k = 1:numel(inj)
    c = remap(W(inj(k)).cells);
    c = c(c>0);
    qi_well(c, k) = sum(state.wellSol(inj(k)).flux(c>0,:), 2);
end

for k = 1:numel(prod)
    c = remap(W(prod(k)).cells);
    c = c(c>0);
    qp_well(c, k) = sum(state.wellSol(prod(k)).flux(c>0,:), 2);
end% distrubute tracer in injector cells according to injector rates
if reverse
    [qi_well, qp_well] = deal(-qi_well, -qp_well);
end
w   = bsxfun(@rdivide, qi_well, sum(qi_well));
tr0 = bsxfun(@rdivide, w, pv(cix));
pv = pv(cix);
end