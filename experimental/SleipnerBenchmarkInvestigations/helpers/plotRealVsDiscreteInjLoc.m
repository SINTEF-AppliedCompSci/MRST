function [ hfig, hax ] = plotRealVsDiscreteInjLoc(Gt, bf, wellXcoord, wellYcoord, wellCoord_x, wellCoord_y)
% TODO - implement this function to simply add the injection locations to
% an existing (currently open) figure.

    figure; hold on

    % actual location
    plot(wellXcoord,wellYcoord,'ok', ...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor','g',...
            'MarkerSize',10)
        
    % simulated location
    plot(wellCoord_x,wellCoord_y,'xk', ...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor','g',...
            'MarkerSize',10)
        
    legend('Actual injection location','Simulation injection location')

%     rectangle('Position',[zoomX1 zoomY1 zoomX2-zoomX1 zoomY2-zoomY1], 'EdgeColor','r', 'LineWidth',3,...
%               'LineStyle','--')

    plotFaces(Gt, bf, 'EdgeColor','k', 'LineWidth',3);
    axis equal %tight
    
    hfig = gcf;
    hax  = gca;

end

