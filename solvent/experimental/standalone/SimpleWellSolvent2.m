classdef SimpleWellSolvent2 < SimpleWell
    
    methods 
    function wellSol = updateConnectionPressureDrop(well, wellSol0, wellSol, model, q_s, bhp, packed, dt, iteration)
        % Update the pressure drop within the well bore, according to a
        % hydrostatic pressure distribution from the bottom-hole to the
        % individual perforations.
        %
        % To avoid dense linear systems, this update only happens at
        % the start of each nonlinear loop.
        if iteration ~= 1
            return
        end
        [p, mob, rho, dissolved, comp, wellvars] = unpackPerforationProperties(packed);
        for i = 1:numel(rho)
            rho{i} = double(rho{i});
            if ~isempty(dissolved)
                for j = 1:numel(dissolved{i})
                    dissolved{i}{j} = double(dissolved{i}{j});
                end
            end
        end
        rho     = cell2mat(rho);
        active = model.getActivePhases();
        numPh = nnz(active);

        rhoS = model.getSurfaceDensities();
        b = phaseDensitiesTobfactor(rho, rhoS, dissolved);
        
        w = well.W;
        if ~isfield(w, 'topo')
            nperf = numel(w.cells);
            w.topo = [(0:(nperf-1))', (1:nperf)'];
        end
        qs = wellSol.cqs; %volumetric in-flux at standard conds
        C = well.wb2in(w);            % mapping wb-flux to in-flux
        wbqs  = abs(C\qs);       % solve to get well-bore fluxes at surface conds
        wbqst = sum(wbqs, 2);   % total wb-flux at std conds
        % if flux is zero - just use compi
        zi = wbqst == 0;
        if any( zi )
            wbqs(zi,:)  = ones(nnz(zi),1)*w.compi;
            wbqst(zi,:) = sum(wbqs(zi,:), 2);
        end
        % Compute mixture at std conds:
        mixs = wbqs ./ (wbqst*ones(1,numPh));
        % compute volume ratio Vr/Vs
        volRat = well.compVolRat(mixs, p, b, model);
        % Mixture density at connection conds (by using static b's)
        rhoMix = (mixs*rhoS')./volRat;
        % rhoMix is now density between neighboring segments given by
        % topo(:,1)-topo(:,2) computed by using conditions in well-cell
        % topo(:,2). This is probably sufficiently accurate.

        % get dz between segment nodes and bh-node1. This is a simple
        % hydrostatic distribution.
        dpt = [0; w.dZ];
        dz  = diff(dpt);
        g   = norm(gravity);
        ddp = g*rhoMix.*dz; % p-diff between connection neighbors
        % well topology assumes we can traverse from top down, but we
        % use a loop just in case of crazy ordering. If this loop does
        % not converge, the solver will throw an error.
        cdp    = nan(size(ddp));
        cdp(1) = ddp(1);
        its = 0; maxIts = 100;
        while and(any(isnan(cdp)), its<maxIts)
            its = its +1;
            % Traverse from top node and down with the pressure
            % differentials found earlier.
            for cnr = 2:numel(cdp)
                cdp(w.topo(cnr,2)) = cdp(w.topo(cnr,1)) + ddp(cnr);
            end
        end
        if its == maxIts
            % If this loop did not converge, something is wrong with
            % either the densities or the well itself. Regardless of
            % reason, we throw an error.
            error(['Problem with topology for well: ', wellSol.name, '. Segments appear not to be connected'])
        end
        wellSol.cdp = cdp;
    end
    end
end

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
