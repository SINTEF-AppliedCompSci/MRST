function fn = plotLimiter(model, varargin)

    opt = struct('variable', 's'      , ...
                 'plot1d'  , false    , ...
                 'view'    , [100, 20], ...
                 'position', []       , ...
                 'n'       , 200      , ...
                 'zlim'    , [0,1]    , ...
                 'pbaspect', [1,1,0.3]  );
    [opt, patchOpt] = merge_options(opt, varargin{:});
    % Get coordinates for plotting
    coords = getPlotCoordinates(model.G, 'n', opt.n, 'plot1d', opt.plot1d);
    % Set figure position
    if opt.plot1d
        opt.position = [0,0,800,400];
    else
        opt.position = [0,0,800,800];
    end
    % Make figure handle
    fig = figure('Position', opt.position);
    disc = model.discretization;
    % Make plot functions
    if opt.plot1d
        [hu, su] = plotSaturationDG(disc, [], 'plot1d', true, 'coords', coords, 'makePlot', false);
        [hl, sl] = plotSaturationDG(disc, [], 'plot1d', true, 'coords', coords, 'makePlot', false);
    else
        [hu, su] = plotSaturationDG(disc, [], 'coords', coords, 'makePlot', false);
        [hl, sl] = plotSaturationDG(disc, [], 'coords', coords, 'makePlot', false);
    end
    fn = @(model, states, reports, solver, schedule, simtime) ...
        afterStepFunction(model, states, reports, solver, schedule, simtime, disc, coords, fig, {hu, hl}, {su, sl}, opt, patchOpt);
end

function [model, states, reports, solver, ok] ...
    = afterStepFunction(model, states, reports, solver, schedule, simtime, disc, coords, fig, h, sat, opt, patchOpt)
    computed = cellfun(@(x) ~isempty(x), states);
    current = find(computed, 1, 'last');
    st      = states{current};
    if ishandle(fig)
        set(0, 'CurrentFigure', fig);
    else
        figure();
%         error('Figure handle deleted');
    end
    if current == 1
        clf;
    end
    ax = fig.Children(isa(fig.Children, 'matlab.graphics.axis.Axes'));
    if opt.plot1d
        sul = sat{1}(st.ul);
        sl  = sat{2}(st);
        if isempty(ax)
            x = coords.points(:,1); 
            hold on
            plot(x(:,1), sul, 'color', 'k', 'lineWidth', 2);
            plot(x(:,1), sl , 'color', 'k', 'lineWidth', 4, 'LineStyle', '--');
            ylim([-0.2, 1.2]);
            hold off
        else
            l = ax.Children(arrayfun(@(c) isa(c, 'matlab.graphics.chart.primitive.Line'), ax.Children));
            l(2).YData = sul;
            l(1).YData = sl;
        end
    else
        clf
        subplot(1,2,1)
        plotSaturationDG(disc, st.ul, patchOpt{:}, 'coords', coords);
        view(opt.view)
        camlight
        lighting gouraud
        axis tight
        zlim(opt.zlim)
        pbaspect(opt.pbaspect);
        subplot(1,2,2)
        plotSaturationDG(disc, st, patchOpt{:}, 'coords', coords);
        view(opt.view)
        camlight
        axis tight
        zlim(opt.zlim)
        pbaspect(opt.pbaspect);
        lighting gouraud
    end
    ok = true;
    drawnow();
end