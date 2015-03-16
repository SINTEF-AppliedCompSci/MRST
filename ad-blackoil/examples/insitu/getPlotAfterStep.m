function fn = getPlotAfterStep(state0, model, schedule, varargin)
    opt = struct('plotWell', true, 'plotReservoir', true);
    opt = merge_options(opt, varargin{:});
    
    G = model.G;
    
    W = schedule.control(1).W;
    
    if ~isempty(W) && opt.plotWell
        ws = initWellSolAD(schedule.control(1).W, model, state0);
        sources = {vertcat(ws.qWs), vertcat(ws.qOs), vertcat(ws.qGs)};
        
        model.wellmodel.W = W;
        ws = model.wellmodel.updateWellSolStatistics(ws, sources, model);

        ws0 = {ws; ws};
        
        [hwell, injectWell] = plotWellSols(ws0);
    else
        [hwell, injectWell] = deal(nan);
    end
    
    if opt.plotReservoir
        hdata = figure;
        [~, injData] = plotToolbar(G, {state0; state0});
        axis tight
    else
        [hdata, injData] = deal(nan);
    end
    hdata = double(hdata);
    hwell = double(hwell);
    fn = @(model, states, reports, schedule, simtime) afterStepFunction(model, states, reports, schedule, simtime, injData, injectWell, hdata, hwell);
end

function [model, states, reports, ok] = afterStepFunction(model, states, reports, schedule, simtime, injData, injectWell, hdata, hwell)
    computed = cellfun(@(x) ~isempty(x), states);
    current = find(computed, 1, 'last');
    
    st = states(computed);   
    rep = reports(computed);
    simtime = simtime(computed);
    if ~isnan(hwell) && ishandle(hdata)
        injData(model.G, st, current);
    end
    
    if ~isnan(hwell) && ishandle(hwell)
        set(0, 'CurrentFigure', hwell);
        ws = cellfun(@(x) x.wellSol, st, 'uniformoutput', false);
        % Note: We are not actually considering the case where ministeps
        % are being inputed, which would result in an error
        T  = cumsum(schedule.step.val(computed));
        if numel(T) ~= numel(ws)
            T = (1:numel(ws))';
        end
        injectWell({ws}, T)
    end
    
    ok = true;
    if 1
        ok = ok & simulationRuntimePanel(model, reports, schedule, simtime);
    end
    

end
