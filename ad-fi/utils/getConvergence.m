function [converged, CNV, MB] = getConvergence(state, eqs, fluid, system, dt)
    % Compute convergence based on total mass balance (tol_mb) and maximum
    % residual mass balance (tol_cnv).

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

    ac = system.activeComponents;
    tol_mb = system.nonlinear.tolMB;
    tol_cnv = system.nonlinear.tolCNV;

    activePhases = [ac.oil, ac.water, ac.gas];

    pv = system.s.pv;
    pvsum = sum(pv);
    %  Kan bruke snittrykk her. Profilering tilsier at getConvergence
    %  koster veldig lite.

    nc = numel(state.pressure);
%     fprintf('%2.5f %2.5f \n', tol_mb, tol_cnv);
    % WATER
    water_ind=sum(activePhases(1:2));
    oil_ind=sum(activePhases(1));
    gas_ind=sum(activePhases);

    if system.activeComponents.water,
        %BW = fluid.BW(state.pressure);
        BW = 1./fluid.bW(state.pressure);
        RW = eqs{water_ind}.val;
        BW_avg = sum(BW)/nc;
        CNVW = BW_avg*dt*max(abs(RW)./pv);
    else
        BW_avg = 0;
        CNVW   = 0;
        RW     = 0;
    end

    % OIL
    if system.activeComponents.oil,
        if system.activeComponents.disgas
            % If we have liveoil, BO is defined not only by pressure, but
            % also by oil solved in gas which must be calculated in cells
            % where gas saturation is > 0.
            BO = 1./fluid.bO(state.pressure, state.rs, state.s(:,gas_ind)>0);
        else
            %BO = fluid.BO(state.pressure);
            BO = 1./fluid.bO(state.pressure);
        end
        RO = eqs{oil_ind}.val;
        BO_avg = sum(BO)/nc;
        CNVO = BO_avg*dt*max(abs(RO)./pv);
    else
        BO_avg = 0;
        CNVO   = 0;
        RO     = 0;
    end

    % GAS
    if system.activeComponents.gas,
        if system.activeComponents.vapoil
            BG = 1./fluid.bG(state.pressure, state.rv, state.s(:,2)>0); % need to fix index...
        else
            BG = 1./fluid.bG(state.pressure);
        end
        RG = eqs{gas_ind}.val;
        BG_avg = sum(BG)/nc;
        CNVG = BG_avg*dt*max(abs(RG)./pv);
    else
        BG_avg = 0;
        CNVG   = 0;
        RG     = 0;
    end

    % Check if material balance for each phase fullfills residual
    % convergence criterion
    MB = abs([BO_avg*sum(RO), BW_avg*sum(RW) BG_avg*sum(RG)]);
    converged_MB  = all(MB < tol_mb*pvsum/dt);

    % Check maximum normalized residuals (maximum mass error)
    CNV = [CNVO CNVW CNVG] ;
%     [CNVO CNVW CNVG]
    converged_CNV = all(CNV < tol_cnv);

%     figure(1); clf; hold on;
%     plot([BO_avg*dt*abs(RO)./pv BW_avg*dt*abs(RW)./pv]);
%     plot(tol_cnv*ones(10,1),'red');
%     plot((1./pv)*min(pv),'green');
%     drawnow;
%     pause(.1)
%     figure(2); plot(pv); drawnow;

    converged = converged_MB & converged_CNV;
%     semilogy(CNV);drawnow
end

% function [CNV BR_avg] = EclipseConvergence(BX, RX, pv, nc, active)
%     if active;
%         BX_avg = sum(BX)/nc;
%         CNV = BX_avg*dt*max(RX./pv);
%     else
%         BR_avg = 0;
%         CNV   = 0;
%     end
% end
