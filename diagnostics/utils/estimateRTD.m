function [dist] = estimateRTD(pv, D, WP, t_end, varargin)
opt = struct('injectorIx', [], ...
             'producerIx', [], ...
             'nbins',      25, ...
             'match_allocation', false);
         
opt = merge_options(opt, varargin{:});

[producerIx, injectorIx] = deal(opt.producerIx, opt.injectorIx);
if ~isempty(producerIx) && ~isempty(injectorIx)
    assert(isfield(D, 'ptof') && isfield(D, 'itof'), ...
            'RTD estimation for subset of all wells require input of individual well TOFs (fields itof/ptof)');
end

if isempty(producerIx), producerIx = (1:numel(D.prod))'; end
if isempty(injectorIx), injectorIx = (1:numel(D.inj))'; end

nreg = numel(producerIx).*numel(injectorIx);

% create output structure
dist = struct('pairIx', zeros(nreg, 2), 't', zeros(opt.nbins+1, nreg), ...
              'volumes', zeros(nreg, 1), 'allocations',  zeros(nreg,1), ...
              'values', zeros(opt.nbins+1, nreg) );              

          
ix = 0;          
for ik = 1:numel(injectorIx)
    for pk = 1:numel(producerIx)
        [pix, iix] = deal(producerIx(pk), injectorIx(ik));
        ix = ix +1;
        % collect data from WP
        dist.pairIx(ix,:) = [iix, pix];
        ixWP = WP.pairIx(:,1)==iix & WP.pairIx(:,2) == pix;
        dist.volumes(ix)      = WP.vols(ixWP);
        dist.allocations(ix)  = sum(WP.inj(iix).alloc(:,pix));
        
        % collect relevant tracer values
        c  = [D.itracer(:, iix), D.ptracer(:, pix)];
        cp = prod(c,2);
        % index to relevant subset
        sub  = cp > 1e-5;
        nsub = nnz(sub);
        
        if nsub > 0
            % compute tof as weighted average of well-tofs
            tof = zeros(nsub, 2);
            if isfield(D, 'itof')
                tof(:,1) = D.itof(sub, iix);
            else
                tof(:,1) = D.tof(sub,1);
            end
            if isfield(D, 'ptof')
                tof(:,2) = D.ptof(sub,pix);
            else
                tof(:, 2) = D.tof(sub,2);
            end
            
            % total tof/residence time
            ttof  = sum(tof, 2);
            % relevant pore volumes
            pvs   = pv(sub).*cp(sub);
            % sort according to ttof
            [ts, order] = sort(ttof);
            % approximate flux through each cell
            flux = pvs(order)./ts;
            
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
                fac   = sum(binflux)/dist.allocations(ix);
                unitbinflux = unitbinflux*fac;
                dispif(mrstVerbose, 'Distribution scaled by %3.2f to match allocation.\n', fac);
            end
            
            % ommit last entry (t_end to infinity)
            dist.t(:, ix)      = edges(i(1:end-1));
            dist.values(:, ix) = unitbinflux(i(1:end-1));
        end
    end
end