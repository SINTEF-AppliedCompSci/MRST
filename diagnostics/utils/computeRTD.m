function dist = computeRTD(state, G, pv, D, WP, W, t_end, varargin)
opt = struct('injectorIx', [], ...
             'producerIx', [], ...
             'nsteps',     50);
opt = merge_options(opt, varargin{:});

[pix, iix] = deal(opt.producerIx, opt.injectorIx);         
if isempty(pix), pix = (1:numel(D.prod))'; end
if isempty(iix), iix = (1:numel(D.inj))'; end       

% index into well-struct
[inj, prod] = deal(D.inj(iix), D.prod(pix));

% extract relevat subset
sub = sum(D.itracer(:, iix), 2) > 1e-5;

dt = t_end/opt.nsteps;
[A, qi_well, qp_well, tr0] = setupSystemComponents(state, G, pv, W, dt, sub, inj, prod);


tr = tr0;
vals = cell(1, numel(iix));
for k = 1:numel(vals)
    vals{k} = zeros(opt.nsteps+1, numel(pix));
end
h = waitbar(0, 'Computing distribution(s) ...');
for k = 1:opt.nsteps
    tr = A\tr;
    for tn = 1:numel(vals)
        curvals  = -(tr(:, tn)'*qp_well);
        curvals(~isfinite(curvals)) = 0;
        vals{tn}(k+1,:) = curvals;
        waitbar(k/opt.nsteps, h);
    end
    %
end
close(h);
t = (0:opt.nsteps)'*(t_end/opt.nsteps);

% create output struct
nreg = numel(pix).*numel(iix);
dist = struct('pairIx', zeros(nreg, 2), 't', zeros(opt.nsteps+1, nreg), ...
              'volumes', zeros(nreg, 1), 'allocations',  zeros(nreg,1), ...
              'values', zeros(opt.nsteps+1, nreg) ); 
ix = 0;          
for ik = 1:numel(iix)
    for pk = 1:numel(pix)
        ix = ix +1;
        % collect data from WP
        dist.pairIx(ix,:) = [iix(ik), pix(pk)];
        ixWP = WP.pairIx(:,1)==iix(ik) & WP.pairIx(:,2) == pix(pk);
        dist.volumes(ix)      = WP.vols(ixWP);
        dist.allocations(ix)  = sum(WP.inj(iix(ik)).alloc(:,pix(pk)));
        dist.t                = repmat(t, [1, nreg]);
        dist.values(:, ix)    = vals{ik}(:, pk);
    end
end
end
    
function [A, qi_well, qp_well, tr0, pv] = setupSystemComponents(state, G, pv, W, dt, sub, inj, prod)
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
A = speye(ncr) - dt*A;

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
tr0 = bsxfun(@rdivide, qi_well, pv(cix));
pv = pv(cix);
end