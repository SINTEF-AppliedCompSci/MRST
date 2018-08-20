function [t, vals] = estimateRTD(pv, D, WP, t_end, varargin)
opt = struct('injectorIx', [], ...
             'producerIx', [], ...
             'nbins',      25, ...
             'match_allocation', true);
opt = merge_options(opt, varargin{:});

[pix, iix] = deal(opt.producerIx, opt.injectorIx);
if ~isempty(pix) && ~isempty(iix)
    assert(isfield(D, 'ptof') && isfield(D, 'itof'), ...
            'RTD estimation for subset of all wells require input of individual well TOFs (fields itof/ptof)');
end

if isempty(pix), pix = (1:numel(D.prod))'; end
if isempty(iix), iix = (1:numel(D.inj))'; end

% collect relevant tracer values
c  = [sum(D.itracer(:, iix), 2), sum(D.ptracer(:, pix), 2)];
cp = prod(c,2);
% index to relevant subset 
sub  = cp > 1e-5;
nsub = nnz(sub);

if nsub == 0
    warning('Inactive well, empty distribution')
     [t, vals] = deal([]);
     return;
end
% compute tof as weighted average of well-tofs
tof = zeros(nsub, 2);
if isfield(D, 'itof')
    itof = D.itof(sub,iix);
    itof(isinf(itof)) = 1e5*year;
    tof(:,1) = sum(D.itracer(sub, iix).*itof ,2)./c(sub, 1);
else
    tof(:,1) = D.tof(sub,1);
end
if isfield(D, 'ptof')
    ptof = D.ptof(sub,pix);
    ptof(isinf(ptof))=1e5*year;
    tof(:,2) = sum(D.ptracer(sub, pix).*ptof ,2)./c(sub, 2);
else
    tof(:, 2) = D.tof(sub,2);
end

% total tof/residence time
ttof  = sum(tof, 2);
% relevant pore volumes
pv   = pv(sub).*cp(sub);
% sort according to ttof
[ts, order] = sort(ttof);
% approximate flux through each cell
flux = pv(order)./ts;

% bin edges
edges = [(0:opt.nbins)*(t_end/opt.nbins), inf]; 
[~, ~, bins] = histcounts(ts, edges);
% sum fluxes for each bin
binflux = accumarray(bins, flux);
% divide by bin-length to get unit flux
i = 1:min(numel(binflux)+1,numel(edges));
unitbinflux = [0; binflux./diff(edges(i)).'];

% normalize so total flux equals allocation
if opt.match_allocation
    alloc = sum(arrayfun(@(x)sum(sum(x.alloc(:,pix))), WP.inj(iix)));
    fac   = sum(binflux)/alloc;
    unitbinflux = unitbinflux*fac;
    dispif(mrstVerbose, 'Distribution scaled by %3.2f to match allocation.\n', fac);
end

% ommit last entry (t_end to infinity)
[t,vals] = deal(edges(i(1:end-1)), unitbinflux(i(1:end-1)));
end