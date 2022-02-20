function [dist] = estimateRTD(pv, D, WP, varargin)
%Estimate residence time distributions based on computed tof-values
%
% SYNOPSIS:
%   dist = estimateRTD(pv, D, WP, 'pn1', pv1, ...)
% 
% DESCRIPTION:
%  This function estmates well-pair RTDs based on computed TOFs and tracer
%  fields. The RDT approximates an imaginary tracer which ditributes equally 
%  among all phases (follows the total flux-field). A tracer pulse is 
%  injected at injector i at time zero and produced at producer p. The RTD
%  is scaled such that
%     RTD_{ip}(t) ~= (rate producer p / total tracer mass) * c_p(t), 
%  where c_p(t) is the tracer concentration in producer p. 
% 
%  With the above scaling the RTD has units [s]^-1
%   * the integral of the RTD approximates fractional recovery 
%    (produced mass/injected mass) 
%   * the mean of the RTD time i-to-p allocation approximates the i-to-p 
%     allocation volume
%
% RETURNS:
%  dist - structure with fields
%    pairIx            nreg x 2 each row index to injector/producer 
%    t                 nbin x 1 discrete times
%    values            nbin x 1 discrete RTD-values
%    volumes           nreg x 1 interaction volume for each well pair
%    allocations       nreg x 1 interaction allocation for each well pair
%    injectorFlux      ninj injector total rates
%    producerFlux      nprod producer total rates
%
% SEE ALSO
%  computeRTD

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
             'nbins',      100, ...
             'match_allocation', true);
         
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
dist = struct('pairIx',            nan(nreg, 2), ...
              't',                 nan(opt.nbins, nreg), ...
              'volumes',           nan(nreg, 1), ...
              'allocations',       nan(nreg,1), ...
              'values',            nan(opt.nbins, nreg), ...
              'injectorFlux',      nan(numel(injectorIx),1), ...
              'producerFlux',      nan(numel(producerIx),1));
% 
for ik = 1:numel(injectorIx)
    dist.injectorFlux(ik) = sum( sum(WP.inj(injectorIx(ik)).alloc) );
end

for pk = 1:numel(producerIx)
    dist.producerFlux(pk) = sum( sum(WP.prod(producerIx(pk)).alloc) );
end          
          
          
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
        
        q_inj = sum(sum(WP.inj(iix).alloc));
        % collect relevant tracer values
        c  = [D.itracer(:, iix), D.ptracer(:, pix)];
        cp = prod(c,2);
        % index to relevant subset
        sub  = cp > 1e-5;
        nsub = nnz(sub);
        
        if nsub==0, continue; end
        
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
        
        % bin edges: histcounts was introduced in MATLAB R2014b, so we
        % provide a backup solution for older versions
        try
            [~,edges, bins] = histcounts(log10(ts), opt.nbins);
        catch
           t = log10(ts);
           edges = linspace(min(t),max(t+0.01), opt.nbins+1);
           [~,bins] = histc(t, edges);
        end
        edges = 10.^edges;
           
        % sum fluxes for each bin
        binflux = accumarray(bins, flux);
        totflux = sum(binflux);
        if numel(binflux)+1<numel(edges)
            binflux(end+1:numel(edges)-1)=nan;
        end
        % divide by bin-length to get unit flux
        unitbinflux = [0; binflux./diff(edges)'];
        
        % normalize so total flux equals allocation
        if opt.match_allocation
            fac   = totflux/dist.allocations(ix);
            unitbinflux = unitbinflux/fac;
            dispif(mrstVerbose, 'Distribution scaled by %3.2f to match allocation.\n', fac);
        end
        
        % ommit last entry
        i = (1:numel(binflux)).';
        dist.t(i, ix)      = edges(i);
        % scale by q_inj such that distribution represents 1kg injected tracer
        dist.values(i, ix) = unitbinflux(i)/q_inj;
    end
    dist.creator = mfilename;
end
