function plotFormationPressureChanges( states, initPressure, varargin )

localopt.figname = [];
localopt.P_lim   = []; % A set/computed overpressure limit (in Pascals), 
                       % to avoid approaching overburden (fracture) pressure
localopt = merge_options( localopt, varargin{:} );

    for i = 1:numel(states)
        maxp(i) = max( states{i}.pressure );
        minp(i) = min( states{i}.pressure );
        avgp(i) = mean( states{i}.pressure );
    end

    maxinitp = max( initPressure );
    mininitp = min( initPressure );
    avginitp = mean( initPressure );

    figure
    if ~isempty(localopt.figname)
        set(gcf,'name',localopt.figname);
    end
    hold on
    plot([0 1:1:numel(states)], [maxinitp maxp]/1e6, 'LineWidth',3)
    plot([0 1:1:numel(states)], [mininitp minp]/1e6, 'LineWidth',3)
    plot([0 1:1:numel(states)], [avginitp avgp]/1e6, 'LineWidth',3)
    plot([0 1:1:numel(states)], repmat(avginitp/1e6,[1 numel(states)+1]), '--', 'LineWidth',3)
    if ~isempty(localopt.P_lim)
        plot([0 1:1:numel(states)], repmat((localopt.P_lim+avginitp)/1e6,[1 numel(states)+1]), '--', 'LineWidth',3)
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