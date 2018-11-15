classdef SummaryViewer < handle
    % Simple summary viewer. Usage (name vithout extension)
    % >> SummaryViewer(readEclipseSummaryUnFmt([filepath, name]))
    properties
        Figure
        Axes
        selector
        smry
        layout = struct('menuWidth', 200)
    end
    
    methods
        function d = SummaryViewer(smry)
            if isstring(smry)
                d.smry = readEclipseSummaryUnFmt(smry);
            else
                d.smry = smry;
            end
            d.Figure = figure;
            d.Axes   = axes('Parent', d.Figure, 'Units', 'pixels');
            d.selector = SummarySelector(d.smry, 'Parent', d.Figure, 'Units', 'pixels');
            d.selector.panel.ButtonDownFcn = '';
            % don't show controls we dont use
            c = d.selector.controls;
            for k1 = 2:numel(c)
                for k2 = 1:numel(c{k1})
                    c{k1}{k2}.Visible = 'off';
                end
            end
            d.Figure.SizeChangedFcn = @d.updateLayout;
            d.selector.Callback = @d.summaryCallback;
            d.updateLayout();
        end
        
        function summaryCallback(d, src, event)
            [nms, prps] = deal(d.selector.curNames, d.selector.curProps);

            ax = d.Axes;
            if ~isempty(prps)
                cla(ax, 'reset');
                leg = {};
                hold(ax, 'on')
                for k = 1:numel(nms)
                    for l = 1:numel(prps)
                        if any(strcmp(prps{l}, d.smry.getKws(nms{k})))
                            plot(ax, d.selector.time, d.smry.get(nms{k}, prps{l}, ':'), 'LineWidth', 2);
                            leg = [leg, {[nms{k},' - ', prps{l},' [', strtrim(d.smry.getUnit(nms{k}, prps{l})), ']']}]; %#ok
                        end
                    end
                    if ~isempty(leg)
                        %xlabel('time [years]')
                        legend(ax, leg, 'Interpreter', 'none')
                        xtickangle(ax, 30)
                    end
                end
            end
        end
        
        function updateLayout(d, src, event)
            fpos = d.Figure.Position;
            mw = d.layout.menuWidth;
            d.selector.Position = [0, -8, mw, fpos(4)-10];
            d.Axes.Position     = [mw + 50, 50, fpos(3)-mw-75, fpos(4)-75];
        end
    end
end