function fn = getPlotAfterStep(state0, model, schedule)
    G = model.G;
    
    W = schedule.control(1).W;
    
    if ~isempty(W)
        ws = initWellSolAD(schedule.control(1).W, model, state0);
        ws0 = {ws; ws};
        [hwell, injectWell] = plotWellSols(ws0);
    else
        [hwell, injectWell] = deal(nan);
    end
    hdata = figure;
    [~, injData] = plotToolbar(G, {state0; state0});
    fn = @(model, states, reports, schedule) afterStepFunction(model, states, reports, schedule, injData, injectWell, hdata, hwell);
end

function [model, states, reports, ok] = afterStepFunction(model, states, reports, schedule, injData, injectWell, hdata, hwell)
    computed = cellfun(@(x) ~isempty(x), states);
    current = find(computed, 1, 'last');
    
    st = states(computed);    
    if ishandle(hdata)
        injData(model.G, st, current);
    end
    
    if ~isnan(hwell) && ishandle(hwell)
        figure(hwell);
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
end
