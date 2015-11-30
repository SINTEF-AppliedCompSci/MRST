function [ hfig ] = plotOptWellRates( Gt, ta, init, optim, rhoG )



    wcinx   = [init.schedule.control(1).W.cells];
    inrates = [init.schedule.control(1).W.val];
    oprates = [optim.schedule.control(1).W.val];

    % location of wells in formation
    plotWells(Gt, ta, wcinx)

    % compare init and optim well rates
    compareInitOptimRates( inrates , oprates, rhoG )
    



end


function plotWells(Gt, ta, wcinx)

    figure; %set(gcf,'Position',[1 1 1242 825]);

    title('Injectors')
    mapPlot(gcf, Gt, 'traps', ta.traps, 'trapalpha', 0.2, ...
        'rivers', ta.cell_lines, 'rivercolor', [1 0 0], ...
        'maplines', 20, 'wellcells', wcinx, 'well_numbering', true);
    
    % adjust plots
    hfig = gcf;
    set(findobj(hfig.Children,'Type','axes'),'Fontsize',16,'box','on')
    axis(findobj(hfig.Children,'Type','axes'),'equal','tight')

end


function compareInitOptimRates( inrates, oprates, rhoG )

    figure;
    
    title('Initial and Optimized Rates')
    
    if numel(inrates)==1
        % hack for bar plotting
        inrates = [inrates, 0];
        oprates = [oprates, 0];
        bar([inrates; oprates]')
        legend('unoptimized','optimized')
        ylabel('Rates (m^3/s)')
        xlabel('Well Number')
        % adjust XLim to focus on one bin
        hfig = gcf;
        set(findobj(hfig.Children,'Type','axes'), 'XLim',[0.5 1.5])
    else
        %bar([inrates; oprates]')
        [AX,H1,H2] = plotyy(1:numel(inrates), [inrates; oprates]', ...
                            1:numel(inrates), [inrates; oprates]'.*rhoG.*(1/1e9).*convertFrom(1,year), ...
                            'bar','bar');

        % ensure yy-axis limits are at same scale
        ax1 = AX(1);
        ax2 = AX(2);
        set(ax2,'YLim',[ax1.YLim(1) ax1.YLim(2)*rhoG*(1/1e9)*convertFrom(1,year)])

        legend('unoptimized','optimized')
        ax1.YLabel.String = 'Rates (m^3/s)';
        ax2.YLabel.String = 'Rates (Mt/yr)';
        xlabel('Well Number')
    end
    
    % adjust plots
    hfig = gcf;
    set(findobj(hfig.Children,'Type','axes'),'Fontsize',16,'box','on')
    
end