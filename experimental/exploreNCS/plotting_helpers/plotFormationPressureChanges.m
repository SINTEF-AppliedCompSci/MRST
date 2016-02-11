function res = plotFormationPressureChanges( states, initPressure, P_over, Gt, schedule, varargin )

localopt.figname = [];
localopt.P_lim   = []; % A set/computed overpressure limit (in Pascals), 
                       % to avoid approaching overburden (fracture) pressure
localopt = merge_options( localopt, varargin{:} );

    %% Pressure changes using max, min, avg over entire formation grid
    
    for i = 1:numel(states)
        maxp(i) = max( states{i}.pressure );
        minp(i) = min( states{i}.pressure );
        avgp(i) = mean( states{i}.pressure );
    end
    
    maxinitp = max( initPressure );
    mininitp = min( initPressure );
    avginitp = mean( initPressure );
    
    %plotPressChanges(maxinitp, mininitp, avginitp, maxp, minp, avgp, localopt);
    
    %% Bar plot of max pressure reached at well cells, compared to their fracture pressures
    wcinx = [schedule.control(1).W.cells];
    tmp = zeros(numel(states), numel(wcinx));
    for i = 1:numel(states)
        tmp(i,:) = states{i}.pressure(wcinx);
    end
    maxP_at_wellCells = max(tmp);
    tmpPlot = [convertTo(P_over(wcinx), mega*Pascal), convertTo(maxP_at_wellCells, mega*Pascal)'];
    tmpPlot = [tmpPlot, convertTo(initPressure(wcinx), mega*Pascal)];
    
    % hack for bar plot, if only 1 well exists
    if numel(wcinx) == 1
        tmpPlot = [tmpPlot; 0 0 0];
    end
    
    figure;
    hbar = bar(tmpPlot);
    set(hbar(1),'DisplayName','Fracture Pressure','FaceColor','r');
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
    figure; set(gcf, 'Position', [1 1 2439 547])
    
    subplot(1,4,1)
    plotCellData(Gt, convertTo(P_over, mega*Pascal), 'EdgeColor','none')
    plotFaces(Gt, boundaryFaces(Gt), 'EdgeColor','k')
    title('Caprock fracture pressure')
    hcb = colorbar; axis equal tight off
    ylabel(hcb, 'MPa', 'fontsize',16, 'rotation',0)
    
    % and then compare the max pressure reached to the fracture pressure
    amount_surpassed = maxP_encountered - P_over;
    amount_surpassed( amount_surpassed < 0 ) = NaN;

    
    subplot(1,4,2);
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
        title({'Max pressure encountered';'(fracture pressure not surpassed)'})
        isFracPressSurpassed = false;
    else
        isFracPressSurpassed = true;
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
        title({'Max pressure encountered';'(red areas surpassed fracture pressure)'})
  
        subplot(1,4,4)
        plotCellData(Gt, convertTo(amount_surpassed, mega*Pascal), 'EdgeColor','none')
        plotFaces(Gt, boundaryFaces(Gt), 'EdgeColor','k')
        title('Amount surpassed fracture pressure')
        hcb = colorbar; axis equal tight off
        ylabel(hcb, 'MPa', 'fontsize',16, 'rotation',0)
    end
    
    
    
    subplot(1,4,3)
    hbar = bar(tmpPlot);
    set(hbar(1),'FaceColor','r', 'EdgeColor','r');
    set(hbar(2),'FaceColor','y', 'EdgeColor','y');
    set(hbar(3),'FaceColor','b', 'EdgeColor','b');
    hl = legend('Fracture Pressure','Max Pressure Encountered','Initial Pressure');
    set(hl, 'Location','NorthOutside')
    xlabel('Well Number')
    ylabel('Pressure (MPa)')
    % hack for x-axes sizing and labeling
    set(gca, 'xlim', [0 numel(wcinx)+1])
    if numel(wcinx) == 1
        set(gca, 'xlim', [0 numel(wcinx)+1], 'XTick', [1:numel(wcinx) []]);
    end
    box
    grid
    
    hfig = gcf;
    set(findall(hfig,'Type','Text'), 'FontSize',16)
    set(findobj(hfig.Children,'Type','Legend'),'FontSize',14)
    set(findobj(hfig.Children,'Type','axes'), 'FontSize',14) %'Ylim',[0 max(maxp/1e6)])

    %% Plot pressure changes over time
    % for cell where fracture pressure was maximally surpassed, or for cell
    % where pressure was closest to its fracture pressure.
    
    % Determine the cell where ...
    if isFracPressSurpassed
        % pressure surpassed P_over the most.
        [val, inx] = max(maxP_encountered - P_over); 
        
        fprintf('The fracture pressure was surpassed by %d MPa (%d bars), in cell %d.\n', ...
            convertTo(val, mega*Pascal), convertTo(val, barsa), inx)
    else
        % pressure was closest to or was equal to its fracture pressure.
        assert( all(P_over - maxP_encountered) >= 0 )
        [val, inx] = min(P_over - maxP_encountered);
        
        if ~isempty(find(wcinx == inx))
           str = ['well ',num2str(find(wcinx == inx))];
        else
           str = ['not a well']; 
        end
        fprintf(['The fracture pressure was not surpassed, but ...\n', ...
            'the pressure of cell %d (%s) was %4.2f percent of its fracture pressure.\n'], ...
            inx, str, (maxP_encountered(inx)/P_over(inx))*100 )
    end
    
    % Plot pressure change over time for that cell:
    for i = 1:numel(states)
        dp_of_cell(i) = states{i}.pressure(inx);
    end
    initp           = initPressure(inx);
    p_over_of_cell  = P_over(inx);
    
    time_yr = convertTo(cumsum(schedule.step.val), year)';
    
    figure;
    hold on
    plot([0 time_yr], [initp dp_of_cell]/1e6, 'LineWidth',3)
    plot([0 time_yr], repmat(p_over_of_cell/1e6,[1 numel(states)+1]), '--', 'LineWidth',3)
    xlim([0 time_yr(end)])
    xlabel('Time (years since start of injection)')
    ylabel('Pressure (MPa)')
    legend(['cell ',num2str(inx),' pressure'],['cell ',num2str(inx),' fracture pressure'])
    hfig = gcf;
    set(findall(hfig,'Type','Text'), 'FontSize',16)
    set(findobj(hfig.Children,'Type','Legend'),'FontSize',14)
    set(findobj(hfig.Children,'Type','axes'), 'FontSize',14) %'Ylim',[0 max(maxp/1e6)])
    box
    grid
    %close


    
end
