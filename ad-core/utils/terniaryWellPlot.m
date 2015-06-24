function terniaryWellPlot(wellSols, T, ix, varargin)
    opt = struct('plotBHP', true);
    opt = merge_options(opt, varargin{:});
    
    ws0 = wellSols{end};
    if nargin < 2
        T = 1:numel(wellSols);
        xname = 'Step #';
    else
        T = T/year;
        xname = 'Time [year]';
    end
    if nargin < 3
        ix = 1:numel(ws0);
    end
    
    if ischar(ix)
        ix = find(strcmpi({ws0.name}, ix));
    end
    
    flds = {'bhp', 'qWs', 'qOs', 'qGs'};
     [data, names] = getWellOutput(wellSols, flds, ix);
    for i = 1:numel(ix)
        figure('Color', 'w');
        hold on
        qw = squeeze(data(:, i, 2:4));
        
        c = 'brg';
        fac = 1;
        for j = 1:3
            if j == 3, fac = 1000; end
            d = abs(qw(:, j))*day;
            plot(T, d/fac, [c(j), '-'], 'linewidth', 2)
        end
        ax = gca;
        set(ax, 'Color','none');
        xlabel(xname);
        legend('Water', 'Oil', 'Gas', 'Location', 'NorthWest');
        ylabel('Rates [m^3/day] (liquid),  [1000m^3/day] (gas)');
        if opt.plotBHP
            axpos = get(ax, 'Position');
            ax2 = axes('Position',axpos,...
                        'XAxisLocation','top',...
                        'YAxisLocation','right',...
                        'Color','none');
            hold on
            plot(T, data(:, i, 1)/mega, 'k-', 'Parent', ax2, 'linewidth', 2);
            set(ax2, 'XTickLabel', []);
            ylabel('Bottom hole pressure [MPa]');
            legend('BHP', 'Location', 'NorthEast');
            axes(ax);
        end
        grid on
        
        
        % Title
        if ws0(ix).sign == -1
            title(['Production data for ', names{i}]);
        else
            title(['Injection data for ', names{i}]);
        end
    end
    
end
