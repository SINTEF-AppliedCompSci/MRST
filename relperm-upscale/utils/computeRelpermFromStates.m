function upscaled_kr = computeRelpermFromStates(states, model_c, model, schedule_c, schedule)
%Compute relative permeability from a fine-scale simulation
%
% SYNOPSIS:
%   upscaled_kr = computeRelpermFromStates(states, model_c, model, schedule_c, schedule)
%
% REQUIRED PARAMETERS:
%   states     - Fine-scale states
%   model_c    - Upscaled coarse model
%   model      - Fine-scale model used to produce states
%   schedule_c - Upscaled schedule with forces corresponding to model_c.
%
%
% RETURNS:
%   upscaled_kr - Struct containing relative permeability tables for
%                 half-faces and well-cell contacts. No processing is done
%                 to ensure monotonicity or finite values.
%
% SEE ALSO:
%   regularizeSaturationFunction, assignUpscaledRelPerm,
%   mergeHalfFaceRelPerm

%{
Copyright 2009-2016 SINTEF ICT, Applied Mathematics.

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

    states = computeProps(model_c, model, states, schedule_c, schedule);
    
    n_ph = model.water + model.oil + model.gas;
    
    
    N = model_c.operators.N;
    n_f = size(N, 1);
    % All wells in MRST schedules are uniform, so this should be ok
    if isfield(schedule_c.control, 'W')
        W = schedule_c.control(1).W;
    else
        W = [];
    end
    n_w = numel(W);
    n_s = numel(states);
    
    upscaled_kr = cell(1, n_ph);
    for i = 1:n_ph
        rps = defaultRelPermStruct(n_s);
        
        r = cell(n_f, 2);
        [r{:}] = deal(rps);
        w = cell(n_w, 1);
        for wNo = 1:n_w
            nperf = numel(W.cells);
            w{wNo} = struct('perf', {cell(nperf, 1)});
            for perfNo = 1:nperf
                w{wNo}.perf{perfNo} = rps;
            end
        end
        upscaled_kr{i}.reservoir = r;
        upscaled_kr{i}.wells = w;
    end
    T = model_c.operators.T;
    for stepNo = 1:n_s
        state = states{stepNo};
        for ph = 1:n_ph
            % Reservoir
            flux = state.iflux(:, ph);
            pot = state.pot(:, ph);
            mu = state.mu(:, ph);
            
            flag = pot <= 0;
            
            muf = model_c.operators.faceUpstr(flag, mu);
            sat = model_c.operators.faceUpstr(flag, state.s(:, ph));
            kr = -muf.*flux./(T.*pot);
            
            % Reduce to meaningful values
            ok = isfinite(kr) & kr >= 0;
            
            faces = find(ok);
            flag = flag(ok);
            kr = kr(ok);
            sat = sat(ok);
            for i = 1:numel(faces)
                sub = 2 - double(flag(i));
                f = faces(i);
                upscaled_kr{ph}.reservoir{f, sub}.S(stepNo) = sat(i);
                upscaled_kr{ph}.reservoir{f, sub}.kr(stepNo) = kr(i);
            end

            % Wells
            if isfield(state, 'wellSol')
                W = schedule_c.control(schedule_c.step.control(stepNo)).W;
                for wno = 1:numel(W)
                    wc = W(wno).cells;
                    WI = W(wno).WI;
                    mu = state.mu(wc, ph);
                    sw = state.s(wc, ph);
                    wf = state.wellSol(wno).flux(:, ph);
                    
                    pw = state.wellSol(wno).cdp + state.wellSol(wno).bhp;
                    dp = (state.pressure(wc) - pw);
                    krw = -mu.*wf./(WI.*dp);
                    
                    ok = isfinite(krw) & krw >= 0;
                    
                    for perfNo = 1:numel(wc)
                        wfi = wf(perfNo);
                        if ~ok(perfNo) || ~(wfi < 0 || (wfi > 0 && W(wno).compi(ph) > 0))
                            continue
                        end
                        % Assume only one perforation in coarse model atm
                        upscaled_kr{ph}.wells{wno}.perf{perfNo}.S(stepNo) = sw(perfNo);
                        upscaled_kr{ph}.wells{wno}.perf{perfNo}.kr(stepNo) = krw(perfNo);
                    end
                    upscaled_kr{ph}.wells{wno}.cells = wc;
                end
            end
        end
    end
    
    for i = 1:n_ph
        for j = 1:n_f
            for k = 1:2
                upscaled_kr{i}.reservoir{j, k} = cleanUpFaceStruct(upscaled_kr{i}.reservoir{j, k});
            end
        end
        for j = 1:n_w
            for k = 1:numel(upscaled_kr{i}.wells{j}.perf)
                upscaled_kr{i}.wells{j}.perf{k} = cleanUpFaceStruct(upscaled_kr{i}.wells{j}.perf{k});
            end
        end
    end
end

function d = cleanUpFaceStruct(d)
    act = isfinite(d.S);
    s = d.S(act);
    kr = d.kr(act);

    [d.S, ix] = sort(s);
    d.kr = kr(ix);
end

function s = defaultRelPermStruct(n_step)
    s = struct('S',  nan(n_step, 1), ...
               'kr', nan(n_step, 1));
end

function states = computeProps(model, model_f, states, schedule, schedule_f)
    fluid = model.fluid;
    gdz = model.getGravityGradient();
    for k = 1:numel(states)
        state = upscaleState(model, model_f, states{k});
        
        
        [pressure, sw, so, sg] = model.getProps(state, 'pressure', 'sw', 'so', 'sg');
        sat = {sw, so, sg};
        act = model.getActivePhases();
        sat = sat(act);
        
        [kr_all, pot, mu, mob, rho, pressures] = deal(cell(1, 3));
        kr = cell(1, numel(sat));
        [kr{:}] = model.evaluateRelPerm(sat);
        
        
        
        index = 1;
        if model.water
            [~, ~, mob{1}, rho{1}, pressures{1}, upcw, pot{1}] = getFluxAndPropsWater_BO(model, pressure, sat{index}, kr{index}, 1, gdz);
            index = index + 1; 
        end
        
        if model.oil
            if isfield(fluid, 'rsSat')
                isSatO = state.s(:,3)>0;
            else
                isSatO = false;
            end
            [~, ~, mob{2}, rho{2}, pressures{2}, upco, pot{2}]  = getFluxAndPropsOil_BO(model, pressure, sat{index}, kr{index}, 1, gdz, state.rs, isSatO);
            index = index + 1; 
        end
        
        if model.gas
            if isfield(fluid, 'rvSat')
                isSatG = state.s(:,2)>0;
            else
                isSatG = false;
            end
            [~, ~, mob{3}, rho{3}, pressures{3}, upcg, pot{3}] = getFluxAndPropsGas_BO(model, pressure, sat{index}, kr{index}, 1, gdz, state.rv, isSatG);
        end
        
        
        for i = 1:numel(kr)
            mu{i} = kr{i}./mob{i};
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
                for ph = size(state.wellSol(i).flux, 2)
                    state.wellSol(wno).flux(:, ph) = accumarray(W(wno).fperf, ws0.flux(:, ph));
                end
                WI = Wf(wno).WI;
                state.wellSol(wno).cdp = accumarray(W(wno).fperf, ws0.cdp.*WI)./accumarray(W(wno).fperf, WI);
            end
        end
        
        states{k} = state;
    end
end
