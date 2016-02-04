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
    
    
    %% Details of when and where max pressure occurred
    
    % time step of when max pressure occurred
    [maxPvalue, timeStep_of_maxP] = max( maxp );
    
    % cell index where max pressure occurred
    [~, cellIndex_of_maxP] = max( states{timeStep_of_maxP}.pressure );
    

    
    % NEW: plot of grid with areas in red that surpass overburden pressure
    figure; set(gcf, 'Position', [1 1 2439 547])
    
%     %plot the fracture pressure of the formation's caprock
%     [ P_over, ~] = computeOverburdenPressure( Gt, rock2D, ...
%                                         seafloor_depth, water_density);
    % passed in as variable
    
    subplot(1,4,1)
    plotCellData(Gt, convertTo(P_over, mega*Pascal), 'EdgeColor','none')
    plotFaces(Gt, boundaryFaces(Gt), 'EdgeColor','k')
    title('Caprock fracture pressure')
    hcb = colorbar; axis equal tight off
    ylabel(hcb, 'MPa', 'fontsize',16, 'rotation',0)
    
    % and then compare the max pressure reached to the fracture pressure
    %amount_surpassed = states{timeStep_of_maxP}.pressure - P_over;
    amount_surpassed = maxP_encountered - P_over;
    amount_surpassed( amount_surpassed < 0 ) = NaN;

    
    subplot(1,4,2)
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
        plotCellData(Gt, convertTo(amount_surpassed, mega*Pascal), find(amount_surpassed>0), ...
            'FaceColor','r', 'FaceAlpha', 0.3, 'EdgeColor','none')
        plotFaces(Gt, boundaryFaces(Gt), 'EdgeColor','k')
        title({'Max pressure encountered';'(red areas surpassed fracture pressure)'})
            
        subplot(1,4,4)
        plotCellData(Gt, convertTo(amount_surpassed, mega*Pascal), 'EdgeColor','none')
        plotFaces(Gt, boundaryFaces(Gt), 'EdgeColor','k')
        title('Amount surpassed fracture pressure')
        hcb = colorbar; axis equal tight off
        ylabel(hcb, 'MPa', 'fontsize',16, 'rotation',0)
    end
    
    % output:
    %res.maxPressure          = max( states{timeStep_of_maxP}.pressure );
    res.maxP_encountered     = maxP_encountered;
    res.maxAmountSurpassed   = max( amount_surpassed );
    res.isFracPressSurpassed = isFracPressSurpassed;
    
    fprintf('The max pressure reached was %d MPa.\n', convertTo(max(maxP_encountered), mega*Pascal))
    if res.isFracPressSurpassed
        fprintf('The fracture pressure was surpassed by %d MPa.\n', convertTo(res.maxAmountSurpassed, mega*Pascal))
    else
        fprintf('The fracture pressure was not surpassed.\n')
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
    
    
    %% Plot pressure changes over time for cell where max pressure occurred
    for i = 1:numel(states)
        dp_of_cell_of_maxP(i) = states{i}.pressure(cellIndex_of_maxP);
    end
    initp                   = initPressure(cellIndex_of_maxP);
    p_over_of_cell_of_maxP  = P_over(cellIndex_of_maxP);
    
    figure;
    hold on
    plot([0 1:1:numel(states)], [initp dp_of_cell_of_maxP]/1e6, 'LineWidth',3)
    plot([0 1:1:numel(states)], repmat(p_over_of_cell_of_maxP/1e6,[1 numel(states)+1]), '--', 'LineWidth',3)
    xlabel('state number (0 is initState)')
    ylabel('pressure (MPa)')
    legend('cell pressure','cell fracture pressure')
    hfig = gcf;
    set(findall(hfig,'Type','Text'), 'FontSize',16)
    set(findobj(hfig.Children,'Type','Legend'),'FontSize',14)
    set(findobj(hfig.Children,'Type','axes'), 'FontSize',14) %'Ylim',[0 max(maxp/1e6)])
    
    box
    grid
    close


    
end


function plotPressChanges(maxinitp, mininitp, avginitp, maxp, minp, avgp, localopt)
% varargins are arrays of size 1 x numel(states)

    nsts = numel(maxp);
    
    figure
    if ~isempty(localopt.figname)
        set(gcf,'name',localopt.figname);
    end
    hold on
    plot([0 1:1:nsts], [maxinitp maxp]/1e6, 'LineWidth',3)
    plot([0 1:1:nsts], [mininitp minp]/1e6, 'LineWidth',3)
    plot([0 1:1:nsts], [avginitp avgp]/1e6, 'LineWidth',3)
    plot([0 1:1:nsts], repmat(avginitp/1e6,[1 nsts+1]), '--', 'LineWidth',3)
    if ~isempty(localopt.P_lim)
        plot([0 1:1:nsts], repmat((localopt.P_lim+avginitp)/1e6,[1 nsts+1]), '--', 'LineWidth',3)
    end
    xlabel('state number (0 is initState)')
    ylabel('pressure (MPa)')
    if ~isempty(localopt.P_lim)
        legend('max fm press.','min fm press.','mean fm press.', ...
        'reference press. (mean)','fm press. limit')
    else
        legend('max fm press.','min fm press.','mean fm press.', ...
        'reference press. (mean)')
    end
    %title({ ['Formation pressure changes using ']; ...
    %        [ opt.bdryType ' boundary condition']})
        
    hfig = gcf;
    set(findall(hfig,'Type','Text'), 'FontSize',16)
    set(findobj(hfig.Children,'Type','Legend'),'FontSize',14, 'Location','SouthEast')
    set(findobj(hfig.Children,'Type','axes'), 'FontSize',14) %'Ylim',[0 max(maxp/1e6)])
    
    box
    grid
    
end