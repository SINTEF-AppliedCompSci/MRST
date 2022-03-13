function T_struct = computeTransmissibilityFromStates(p, states, model, schedule, varargin)
%Undocumented Utility Function

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

    opt = struct('reduction', [], ...
                 'relPotentialTol', sqrt(eps), ...
                 'fracFlowTol', 0.02);
    [opt, extra] = merge_options(opt, varargin{:});
         
    model_c = upscaleModelTPFA(model, p);
    schedule_c = upscaleSchedule(model_c, schedule);
    states_c = computeCoarseProps(model_c, model, states, schedule_c, schedule, extra{:});
    
    % compute potential tollerance relative to max reservoir pressure drop
    maxResDp = max(cellfun(@(x)max(vertcat(x.wellSol.bhp))-min(vertcat(x.wellSol.bhp)), states));
    potentialTol = opt.relPotentialTol*maxResDp;
    
    N = model_c.operators.N;
    if isfield(schedule_c.control, 'W')
        W = schedule_c.control(1).W;
        wc = vertcat(W.cells);
    else
        wc = [];
        W = [];
    end
    nw = numel(W);

    nwc = numel(wc);
    nf = size(N, 1);
    n_step = numel(schedule.step.val);
    n_ph = nnz(model.getActivePhases());

    
    WI = zeros(nwc, n_step);
    T = zeros(nf, n_step);
    
    [T_ph, WI_ph] = deal(cell(1, n_ph));
    [T_ph{:}] = deal(T);
    [WI_ph{:}] = deal(WI);
    p2w = getPerforationToWellMapping(W);
    
    for stepNo = 1:n_step
        if isfield(schedule.control, 'W')
            W = schedule.control(schedule.step.control(stepNo)).W;
        end
        state = states_c{stepNo};
        
        % Cell values
        [v, mob, pot] = deal(0);
        mobT = sum(state.mob, 2);
        for ph = 1:n_ph
            v_ph = state.iflux(:, ph);
            mob_ph = state.mob(:, ph);
            pot_ph = state.pot(:, ph);
            
            upc  = pot_ph <= 0;
            mob_ph = model_c.operators.faceUpstr(upc, mob_ph);
            
            T_ph{ph}(:, stepNo) = computeTrans(v_ph, mob_ph, pot_ph, ...
                                  or(abs(pot_ph) < potentialTol, mob_ph./model_c.operators.faceUpstr(upc, mobT) < opt.fracFlowTol) );
            
            v = v + v_ph;
            mob = mob + mob_ph;
            pot = pot + pot_ph;
        end
        pot = pot/n_ph;
        T(:, stepNo) = computeTrans(v, mob, pot, abs(pot) < potentialTol);
        
        % Well values
        if nw > 0
            ws = state.wellSol;
            
            pot = vertcat(ws.pot);
            v = vertcat(ws.flux);
            mob = vertcat(ws.mob);
            
            vT = sum(v, 2);
            mobT = sum(mob, 2);
            
            compi = vertcat(W.compi);
            
            for ph = 1:n_ph
                ci = compi(p2w, ph);
                
                v_ph = v(:, ph);
                mob_ph = mob(:, ph);
                
                isInj = vT > 0;
                
                mob_ph(isInj) = mobT(isInj).*ci(isInj);
                
                WI_ph{ph}(:, stepNo) = computeTrans(v_ph, mob_ph, pot,...
                                  or(abs(pot) < potentialTol, mob_ph./mobT < opt.fracFlowTol) );
            end
            WI(:, stepNo) = computeTrans(vT, mobT, pot, abs(pot) < potentialTol);
        end
    end
    
    W_struct = struct('WI',     WI, ...
                      'WI_ph',  {WI_ph}, ...
                      'wc',     wc, ...
                      'perf2well', p2w);
    R_struct = struct('T_ph',   {T_ph},...
                      'T',      T, ...
                      'N',      N);
    T_struct = struct('reservoir',  R_struct, ...
                      'wells',      W_struct);
end

function t = computeTrans(v, mob, pot, nanflag)
    t = -v./(mob.*pot);
    t(nanflag) = nan;
end
