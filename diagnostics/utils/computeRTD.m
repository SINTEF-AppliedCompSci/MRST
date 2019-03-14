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
opt = struct('injectorIx', [], ...
             'producerIx', [], ...
             'nsteps',     50, ...
             'nbase',       5);
opt = merge_options(opt, varargin{:});

[pix, iix] = deal(opt.producerIx, opt.injectorIx);         
if isempty(pix), pix = (1:numel(D.prod))'; end
if isempty(iix), iix = (1:numel(D.inj))'; end       


% create output struct
nreg = numel(pix).*numel(iix);
dist = struct('pairIx',         nan(nreg, 2), ...
              't',              nan(2*opt.nsteps+1, nreg), ...
              'volumes',        nan(nreg, 1), ...
              'allocations',    nan(nreg,1), ...
              'values',         nan(2*opt.nsteps+1, nreg), ...
              'injectorFlux',   nan(numel(iix),1), ...
              'producerFlux',   nan(numel(pix),1)); 
          

% 
for ik = 1:numel(iix)
    dist.injectorFlux(ik) = sum( sum(WP.inj(iix(ik)).alloc) );
end

for pk = 1:numel(pix)
    dist.producerFlux(pk) = sum( sum(WP.prod(pix(pk)).alloc) );
end

ix = 0;          
for ik = 1:numel(iix)
    for pk = 1:numel(pix)
        ix = ix +1;
        % collect data from WP
        dist.pairIx(ix,:) = [iix(ik), pix(pk)];
        ixWP = WP.pairIx(:,1)==iix(ik) & WP.pairIx(:,2) == pix(pk);
        dist.volumes(ix)      = WP.vols(ixWP);
        dist.allocations(ix)  = sum(WP.inj(iix(ik)).alloc(:,pix(pk)));
    end
end


% index into well-struct
[inj, prod] = deal(D.inj(iix), D.prod(pix));

% extract relevat subset
sub = sum(D.itracer(:, iix), 2) > 1e-5;

% devide into two or more simulation periods
% 1. 0  <= t <   t1 = min(vol/alloc)  , i.e., minimal interaction mean RT
% 2. t1 <= t < nbase*t1
%     |
% n. nbase^(n-2)*t1 <= t <=  maxtof

% only check allocations > total*10^-3
checkIx = dist.allocations > 1e-3*sum(dist.allocations);
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
t2 = min(maxTOF, 100*pvi);

% find number of simulation periods
n = ceil(log(t2/t1)/log(opt.nbase));
t = t2./(opt.nbase.^(n:-1:0));
dt = diff([0 t], 1, 2)./opt.nsteps;

% allocate 
[dist.t, dist.values] = deal(nan((n+1)*opt.nsteps+1, nreg));


%t  = [1*pvi, min(maxTOF, 100*pvi)];
%dt = t/opt.nsteps;

% setup system
[sysmat, qp_well, tr0] = setupSystemComponents(state, G, pv, W, sub, inj, prod);

tr = tr0;
vals = cell(1, numel(iix));
for k = 1:numel(vals)
    vals{k} = zeros((n+1)*opt.nsteps+1, numel(pix));
end
h = waitbar(0, 'Computing distribution(s) ...');
cnt = 0;
for ti = 1:numel(dt)
    A = sysmat(dt(ti));
    for k = 1:opt.nsteps
        cnt = cnt+1;
        tr = A\tr;
        for tn = 1:numel(vals)
            curvals  = -(tr(:, tn)'*qp_well);
            curvals(~isfinite(curvals)) = 0;
            vals{tn}(cnt+1,:) = curvals;
            waitbar(cnt/((n+1)*opt.nsteps), h);
        end
    end
end
close(h);
t = arrayfun(@(x)repmat(x, [opt.nsteps, 1]), dt, 'UniformOutput', false);
t = cumsum([0; vertcat(t{:})]);
dist.t = repmat(t, [1, nreg]);
ix = 0;          
for ik = 1:numel(iix)
    for pk = 1:numel(pix)
        ix = ix +1;
        dist.values(:, ix)    = vals{ik}(:, pk);
    end
end
dist.creator = mfilename;
end
    
function [sysmat, qp_well, tr0, pv] = setupSystemComponents(state, G, pv, W, sub, inj, prod)
nc = G.cells.num;
N  = G.faces.neighbors;
if size(N, 1) ~= numel(state.flux) % flux on interior faces
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
v   = sum(state.flux(fix,:),2);
neg = v < 0;
N(neg,:) = N(neg, [2 1]);
v(neg)   = -v(neg);

% setup production sources
q       = sum(vertcat(state.wellSol.flux), 2);
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
w   = bsxfun(@rdivide, qi_well, sum(qi_well));
tr0 = bsxfun(@rdivide, w, pv(cix));
pv = pv(cix);
end