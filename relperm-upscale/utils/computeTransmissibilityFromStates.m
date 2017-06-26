function T_struct = computeTransmissibilityFromStates(p, states, model, schedule, varargin)
    opt = struct('reduction', [] ...
                 );
    [opt, extra] = merge_options(opt, varargin{:});
    
    model_c = upscaleModelTPFA(model, p);
    schedule_c = upscaleSchedule(model_c, schedule);
    states_c = computeCoarseProps(model_c, model, states, schedule_c, schedule, extra{:});
    
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
        
        for ph = 1:n_ph
            v_ph = state.iflux(:, ph);
            mob_ph = state.mob(:, ph);
            pot_ph = state.pot(:, ph);
            
            upc  = pot_ph <= 0;
            mob_ph = model_c.operators.faceUpstr(upc, mob_ph);
            
            T_ph{ph}(:, stepNo) = computeTrans(v_ph, mob_ph, pot_ph);
            
            v = v + v_ph;
            mob = mob + mob_ph;
            pot = pot + pot_ph;
        end
        pot = pot/n_ph;
        T(:, stepNo) = computeTrans(v, mob, pot);
        
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
                
                WI_ph{ph}(:, stepNo) = computeTrans(v_ph, mob_ph, pot);
            end
            WI(:, stepNo) = computeTrans(vT, mobT, pot);
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

function t = computeTrans(v, mob, pot)
    t = -v./(mob.*pot);
end