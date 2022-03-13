function upscaled_kr = computeRelpermFromStates(states, model_c, model, schedule_c, schedule, varargin)
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
%   `regularizeSaturationFunction`, `assignUpscaledRelPerm`,
%   `mergeHalfFaceRelPerm`

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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
    opt = struct('relPotentialTol', sqrt(eps));
    [opt, extra] = merge_options(opt, varargin{:});
    
    % compute potential tollerance relative to max reservoir pressure drop
    maxResDp = max(cellfun(@(x)max(vertcat(x.wellSol.bhp))-min(vertcat(x.wellSol.bhp)), states));
    potentialTol = opt.relPotentialTol*maxResDp;
    
    states = computeCoarseProps(model_c, model, states, schedule_c, schedule, extra{:});
    
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
            ok = isfinite(kr) & kr >= 0 & abs(pot) > potentialTol;
            
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
                    dp = state.wellSol(wno).pot;
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

