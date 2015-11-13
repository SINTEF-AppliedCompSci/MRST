function [ hfig ] = plotOptWellRates( Gt, ta, init, optim )



    wcinx   = [init.schedule.control(1).W.cells];
    inrates = [init.schedule.control(1).W.val];
    oprates = [optim.schedule.control(1).W.val];

    % location of wells in formation
    plotWells(Gt, ta, wcinx)

    % compare init and optim well rates
    compareInitOptimRates( inrates , oprates )
    



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


function compareInitOptimRates( inrates, oprates )

    figure;
    
    title('Initial and Optimized Rates')
    
    bar([inrates; oprates]')
    
    legend('unoptimized','optimized')
    
    ylabel('Rates (m^3/s)') % or convert to Mt/year !!!
    xlabel('Well Number')

    % adjust plots
    hfig = gcf;
    set(findobj(hfig.Children,'Type','axes'),'Fontsize',16,'box','on')
    
end