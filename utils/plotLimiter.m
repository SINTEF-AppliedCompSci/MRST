function fn = plotLimiter(model, varargin)

    opt = struct('variable', 's'      , ...
                 'plot1d'  , false    , ...
                 'view'    , [100, 20]);
    [opt, patchOpt] = merge_options(opt, varargin{:});
    coords = getPlotCoordinates(model.G, 'n', 200, 'plot1d', opt.plot1d);
    h = figure();
    
    disc = model.discretization;
    fn = @(model, states, reports, solver, schedule, simtime) afterStepFunction(model, states, reports, solver, schedule, simtime, disc, coords, h, opt, patchOpt);
end

function [model, states, reports, solver, ok] = afterStepFunction(model, states, reports, solver, schedule, simtim, disc, coords, h, opt, patchOpt)
    computed = cellfun(@(x) ~isempty(x), states);
    current = find(computed, 1, 'last');
    st      = states{current};
    if ishandle(h)
        set(0, 'CurrentFigure', h);
    else
        figure(h);
    end
    clf;
    if opt.plot1d
        hold on
        plotSaturationDG(disc, st.ul, 'n', 500, 'plot1d', true, 'color', 'k', 'linew', 2);
        plotSaturationDG(disc, st, 'n', 500, 'plot1d', true, 'color', 'k', 'linew', 4, 'LineStyle', '--'); hold off
        ylim([-0.2,1.2]);
    else
        subplot(1,2,1)
        plotSaturationDG(disc, st.ul, patchOpt{:});
        view(opt.view)
        camlight
        lighting gouraud
        zlim([-0.2, 1.2])
        subplot(1,2,2)
        plotSaturationDG(disc, st, patchOpt{:});
        view(opt.view)
        camlight
        zlim([-0.2, 1.2])
        lighting gouraud
    end
    ok = true;
    drawnow();
end