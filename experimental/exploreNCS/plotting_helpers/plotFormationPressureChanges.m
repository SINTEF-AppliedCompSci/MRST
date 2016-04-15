function res = plotFormationPressureChanges( states, initPressure, P_over, plim, Gt, schedule, varargin )

    opt.makePlot = true;
    opt.figname = [];
    opt.outputOn = true;
    opt = merge_options( opt, varargin{:} );
    
    gravity on;

if opt.makePlot
    
    %% Bar plot of max pressure reached at well cells, compared to their fracture pressures
    wcinx = [schedule.control(1).W.cells];
    tmp = zeros(numel(states), numel(wcinx));
    for i = 1:numel(states)
        tmp(i,:) = states{i}.pressure(wcinx);
    end
    maxP_at_wellCells = max(tmp);
    tmpPlot = [convertTo(plim(wcinx), mega*Pascal), convertTo(maxP_at_wellCells, mega*Pascal)'];
    tmpPlot = [tmpPlot, convertTo(initPressure(wcinx), mega*Pascal)];
    
    % hack for bar plot, if only 1 well exists
    if numel(wcinx) == 1
        tmpPlot = [tmpPlot; 0 0 0];
    end
    
    figure;
    hbar = bar(tmpPlot);
    set(hbar(1),'DisplayName','Pressure limit','FaceColor','r');
    set(hbar(2),'DisplayName','Max Pressure Encountered','FaceColor','y');
    set(hbar(3), 'FaceColor','b');
    hl = legend('Fracture Pressure','Max Pressure Encountered','Initial Pressure');
    set(hl, 'Location','NorthOutside')
    xlabel('Well Number')
    ylabel('Pressure (MPa)')
    % hack for x-axes sizing and labeling
    set(gca, 'xlim', [0 numel(wcinx)+1], 'XTick', [1:numel(wcinx) []]);

    hfig = gcf;
    set(findall(hfig,'Type','Text'), 'FontSize',16)
    set(findobj(hfig.Children,'Type','Legend'),'FontSize',14)
    set(findobj(hfig.Children,'Type','axes'), 'FontSize',14) %'Ylim',[0 max(maxp/1e6)])
    
    box
    grid
    close
    
    %% Plot of max pressure encountered during simulation
    tmp = []; for i=1:numel(states), tmp = [tmp, states{i}.pressure]; end
    maxP_encountered = max(tmp')';
    
    
    %% Plot of grid with areas in red that surpass overburden pressure
    %figure; set(gcf, 'Position', [2562 2 2439 547])
    
    %subplot(1,4,1)
    figure; set(gcf,'Position',[2727 42 568 423])
    plotCellData(Gt, convertTo(P_over, mega*Pascal), 'EdgeColor','none')
    plotFaces(Gt, boundaryFaces(Gt), 'EdgeColor','k')
    title('Caprock fracture pressure')
    hcb = colorbar; axis equal tight off
    ylabel(hcb, 'MPa', 'fontsize',16, 'rotation',0)
    hfig = gcf;
    set(findall(hfig,'Type','Text'), 'FontSize',16)
    set(findobj(hfig.Children,'Type','Legend'),'FontSize',14)
    set(findobj(hfig.Children,'Type','axes'), 'FontSize',14) %'Ylim',[0 max(maxp/1e6)])

    
    % and then compare the max pressure reached to the overburden pressure
    amount_surpassed = maxP_encountered - P_over;
    amount_surpassed( amount_surpassed < 0 ) = NaN;

    
    %subplot(1,4,2);
    figure; set(gcf,'Position',[3303 42 568 423])
    plotCellData(Gt, convertTo(maxP_encountered, mega*Pascal), 'EdgeColor','none')
    plotFaces(Gt, boundaryFaces(Gt), 'EdgeColor','k')
    hold on
    xpos = Gt.cells.centroids([schedule.control(1).W.cells],1);
    ypos = Gt.cells.centroids([schedule.control(1).W.cells],2);
    plot(xpos, ypos, 'ok', 'LineWidth',3, 'MarkerFaceColor','k')
    labels = [repmat(' ', numel(xpos), 1), num2str([1:numel(xpos)]')];
    text(xpos, ypos, labels, 'fontsize', 24);
    hcb = colorbar; axis equal tight off
    ylabel(hcb, 'MPa', 'fontsize',16, 'rotation',0)
    if isempty( find(amount_surpassed>0) )
        title({'Max pressure encountered';'(overburden pressure not surpassed)'})
        isOverBurdenPressSurpassed = false;
        hfig = gcf;
        set(findall(hfig,'Type','Text'), 'FontSize',16)
        set(findobj(hfig.Children,'Type','Legend'),'FontSize',14)
        set(findobj(hfig.Children,'Type','axes'), 'FontSize',14) %'Ylim',[0 max(maxp/1e6)])
    else
        isOverBurdenPressSurpassed = true;
        % when plotting cell data, patch error occurs if only one cell data
        % point is to be plotted, thus we pass in the cell data as doubled.
        inx2plot = find(amount_surpassed>0);
        if numel(inx2plot) == 1
            plotCellData(Gt, convertTo(amount_surpassed, mega*Pascal), [inx2plot, inx2plot], ...
            'FaceColor','r', 'FaceAlpha', 0.3, 'EdgeColor','none')
        else
            plotCellData(Gt, convertTo(amount_surpassed, mega*Pascal), inx2plot, ...
            'FaceColor','r', 'FaceAlpha', 0.3, 'EdgeColor','none') 
        end
        plotFaces(Gt, boundaryFaces(Gt), 'EdgeColor','k')
        title({'Max pressure encountered';'(red areas surpassed overburden pressure)'})
        hfig = gcf;
        set(findall(hfig,'Type','Text'), 'FontSize',16)
        set(findobj(hfig.Children,'Type','Legend'),'FontSize',14)
        set(findobj(hfig.Children,'Type','axes'), 'FontSize',14) %'Ylim',[0 max(maxp/1e6)])
        
        %subplot(1,4,4)
        figure; set(gcf,'Position',[4447 39 568 423])
        plotCellData(Gt, convertTo(amount_surpassed, mega*Pascal), 'EdgeColor','none')
        plotFaces(Gt, boundaryFaces(Gt), 'EdgeColor','k')
        title('Amount surpassed overburden pressure')
        hcb = colorbar; axis equal tight off
        ylabel(hcb, 'MPa', 'fontsize',16, 'rotation',0)
        hfig = gcf;
        set(findall(hfig,'Type','Text'), 'FontSize',16)
        set(findobj(hfig.Children,'Type','Legend'),'FontSize',14)
        set(findobj(hfig.Children,'Type','axes'), 'FontSize',14) %'Ylim',[0 max(maxp/1e6)])
    end

    
    
    
    %subplot(1,4,3)
    figure; set(gcf,'Position',[3875 40 568 423])
    hbar = bar(tmpPlot);
    set(hbar(1),'FaceColor','r', 'EdgeColor','r');
    set(hbar(2),'FaceColor','y', 'EdgeColor','k'); % outline in black to help it stand out
    set(hbar(3),'FaceColor','b', 'EdgeColor','b');
    hl = legend('Limit','Max reached','Initial');
    set(hl, 'Location','NorthEast')
    xlabel('Well Number')
    ylabel('Pressure (MPa)')
    % hack for x-axes sizing and labeling
    set(gca, 'xlim', [0 numel(wcinx)+1])
    if numel(wcinx) == 1
        set(gca, 'xlim', [0 numel(wcinx)+1], 'XTick', [1:numel(wcinx) []]);
    end
    grid
    box on;
    hfig = gcf;
    set(findall(hfig,'Type','Text'), 'FontSize',16)
    set(findobj(hfig.Children,'Type','Legend'),'FontSize',14)
    set(findobj(hfig.Children,'Type','axes'), 'FontSize',14) %'Ylim',[0 max(maxp/1e6)])
    
end

    %% Plot pressure changes over time
    % for cell where the maximum percentage of plim was reached
    [perc_of_plim_reach, perc_of_Pover_reach, inx] = ...
                report_maxPercentage_plim_reached( states, plim, P_over, ...
                'outputOn', opt.outputOn);
    res.worst_percent_of_OBpress_reached = perc_of_Pover_reach;
    if opt.outputOn
        if ~isempty(find([schedule.control(1).W.cells] == inx, 1))
            % inx is a well cell
            well_num = find([schedule.control(1).W.cells] == inx, 1);
            fprintf('Cell %d is well %d.\n', inx, well_num)
        else
            % inx is not a well cell 
            fprintf('Cell %d is not a well cell.\n', inx)
        end
    end

if opt.makePlot
    % Plot pressure change over time for that cell:
    plot_p_evol_of_cell(Gt, states, initPressure, P_over, plim, inx, schedule);
end
    
end

function plot_p_evol_of_cell(Gt, states, initPressure, P_over, plim, inx, schedule)

    % Plot pressure change over time for that cell:
    for i = 1:numel(states)
        dp_of_cell(i) = states{i}.pressure(inx);
    end
    initp           = initPressure(inx);
    p_over_of_cell  = P_over(inx);
    p_lim_of_cell   = plim(inx);
    
    time_yr = convertTo(cumsum(schedule.step.val), year)';
    
    figure; set(gcf,'Position',[3126 940 560 420])
    %axes('Position',[0.11 0.15 0.8 0.6])
    hold on
    plot([0 time_yr], [initp dp_of_cell]/1e6, 'LineWidth',3)
    plot([0 time_yr], repmat(p_over_of_cell/1e6,[1 numel(states)+1]), '--', 'LineWidth',3)
    plot([0 time_yr], repmat(p_lim_of_cell/1e6,[1 numel(states)+1]), '--', 'LineWidth',2)
    xlim([0 time_yr(end)])
    xlabel('Time (years since start of injection)')
    ylabel('Pressure (MPa)')
    %legend(['cell ',num2str(inx),' pressure'],['cell ',num2str(inx),' overburden pressure'], ...
    %    ['cell ',num2str(inx),' pressure limit'])
    %legend('P  ','P_o_b  ','P_l_i_m')
    hfig = gcf;
    set(findall(hfig,'Type','Text'), 'FontSize',16)
    %set(findobj(hfig.Children,'Type','Legend'),'FontSize',14,...
    %    'Location','East', 'Orientation','vertical')
    set(findobj(hfig.Children,'Type','axes'), 'FontSize',14) %'Ylim',[0 max(maxp/1e6)])
    box
    grid
    %close
    % Add inset of formation with highlighted cell where pressure reached
    % max percentage of plim
    axes('Position',[0.5 .15 .5 .5]) % start x & y, then axis width & height
    plotGrid(Gt, 'facecolor','none','edgealpha',0.1)
    plotCellData(Gt, ones(Gt.cells.num,1), [inx, inx], 'facecolor','r')
    hold on
    plot(Gt.cells.centroids(inx,1), Gt.cells.centroids(inx,2), 'sq', ...
        'MarkerSize',10, 'LineWidth',2)
    axis equal tight off
    
end
