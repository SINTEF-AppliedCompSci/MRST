function plotFormationPressureChanges( states, var, opt, varargin )

localopt.figname = [];
localopt = merge_options( localopt, varargin{:} );

    for i = 1:numel(states)
        maxp(i) = max( states{i}.pressure );
        minp(i) = min( states{i}.pressure );
        avgp(i) = mean( states{i}.pressure );
    end

    maxinitp = max( var.initState.pressure );
    mininitp = min( var.initState.pressure );
    avginitp = var.ref_p; % mean( var.initState.pressure );

    figure
    if ~isempty(localopt.figname)
        set(gcf,'name',localopt.figname);
    end
    hold on
    plot([0 1:1:numel(states)], [maxinitp maxp]/1e6, 'LineWidth',3)
    plot([0 1:1:numel(states)], [mininitp minp]/1e6, 'LineWidth',3)
    plot([0 1:1:numel(states)], [avginitp avgp]/1e6, 'LineWidth',3)
    plot([0 1:1:numel(states)], repmat(avginitp/1e6,[1 numel(states)+1]), '--', 'LineWidth',3)
    xlabel('state number (0 is initState)')
    ylabel('pressure (MPa)')
    legend('max fm press.','min fm press.','mean fm press.', ...
        'reference press. (mean)')
    title({ ['Formation pressure changes using ']; ...
            [ opt.bdryType ' boundary condition']})
        
    hfig = gcf;
    set(findall(hfig,'Type','Text'), 'FontSize',16)
    set(findobj(hfig.Children,'Type','Legend'),'FontSize',14, 'Location','SouthEast')
    set(findobj(hfig.Children,'Type','axes'), 'FontSize',14,'Ylim',[0 max(maxp/1e6)])
    
    box
    grid
end