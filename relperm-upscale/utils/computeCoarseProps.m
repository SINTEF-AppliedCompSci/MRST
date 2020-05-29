function states = computeCoarseProps(model, model_f, states, schedule, schedule_f, varargin)
%Undocumented Utility Function

%{
Copyright 2009-2020 SINTEF Digital, Mathematics & Cybernetics.

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

    opt = struct('evalmob', true);
    opt = merge_options(opt, varargin{:});
    
    fluid = model.fluid;
    gdz = model.getGravityGradient();
    model = model.validateModel();
    model_f = model_f.validateModel();
    for k = 1:numel(states)
        state = upscaleState(model, model_f, states{k});
        
        
        [pressure, sw, so, sg] = model.getProps(state, 'pressure', 'sw', 'so', 'sg');
        sat = {sw, so, sg};
        act = model.getActivePhases();
        sat = sat(act);
        
        [kr_all, pot, mu, mob, rho, pressures, mu] = deal(cell(1, 3));
        kr = cell(1, numel(sat));
        [kr{:}] = model.evaluateRelPerm(sat);
        
        
        
        index = 1;
        if model.water
            [~, ~, mob{1}, rho{1}, pressures{1}, upcw, pot{1}, mu{1}] = getFluxAndPropsWater_BO(model, pressure, sat{index}, kr{index}, 1, gdz);
            index = index + 1; 
        end
        
        if model.oil
            if isfield(fluid, 'rsSat')
                isSatO = state.s(:,3)>0;
            else
                isSatO = false;
            end
            [~, ~, mob{2}, rho{2}, pressures{2}, upco, pot{2}, mu{2}]  = getFluxAndPropsOil_BO(model, pressure, sat{index}, kr{index}, 1, gdz, state.rs, isSatO);
            index = index + 1; 
        end
        
        if model.gas
            if isfield(fluid, 'rvSat')
                isSatG = state.s(:,2)>0;
            else
                isSatG = false;
            end
            [~, ~, mob{3}, rho{3}, pressures{3}, upcg, pot{3}, mu{3}] = getFluxAndPropsGas_BO(model, pressure, sat{index}, kr{index}, 1, gdz, state.rv, isSatG);
        end
        
        
        if ~opt.evalmob
            pv = model_f.operators.pv;
            p = model.G.partition;
            ix = 1;
            for i = 1:numel(mob)
                if ~isempty(mob{i})
                    mob{i} = accumarray(p, pv.*state.mob(:, ix))./accumarray(p, pv);
                    ix = ix + 1;
                end
            end
        end

        %for i = 1:numel(kr)
        %    mu{i} = kr{i}./mob{i};
        %end
        
        flds = {'kr', 'pot', 'mu', 'mob', 'rho'};
        for i = 1:numel(flds)
            if isfield(state, flds{i})
                state = rmfield(state, flds{i});
            end
        end
        
        [kr_all{act}] = kr{:};
        state = model.setPhaseData(state, kr_all, 'kr');
        state = model.setPhaseData(state, pot, 'pot');
        state = model.setPhaseData(state, mu, 'mu');
        state = model.setPhaseData(state, mob, 'mob');
        state = model.setPhaseData(state, rho, 'rho');
        state = model.setPhaseData(state, pressures, 'pressures');
        state.iflux = state.flux(model.operators.internalConn,:);
        
        if isfield(state, 'wellSol')
            % Assume that well was upscaled using upscaleTrans...
            W = schedule.control(schedule.step.control(k)).W;
            Wf = schedule_f.control(schedule_f.step.control(k)).W;
            for wno = 1:numel(W)
                ws0 = state.wellSol(wno);
                state.wellSol(wno).flux = [];
                for ph = 1:size(ws0.flux, 2)
                    state.wellSol(wno).flux(:, ph) = accumarray(W(wno).fperf, ws0.flux(:, ph));
                end
                WI = Wf(wno).WI;
                wc = W(wno).cells;
                state.wellSol(wno).cdp = accumarray(W(wno).fperf, ws0.cdp.*WI)./accumarray(W(wno).fperf, WI);
                pw = state.wellSol(wno).cdp + state.wellSol(wno).bhp;
                state.wellSol(wno).pot = (state.pressure(wc) - pw);
                % treat injector mobilities later ...
                state.wellSol(wno).mob = state.mob(wc, :);
            end
        end
        
        states{k} = state;
    end
