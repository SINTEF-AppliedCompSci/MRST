function PlotData(data)
%
% DESCRIPTION: plot used for saturation functions data
%
% SYNOPSIS:
%   PlotData(data)
%
% PARAMETERS:
%   data - struct with the plotting data and configurations
%
% RETURNS:
%   plot of the data struct
%
% ----------------------------------
% (c) 2020-2022
% Siroos Azizmohammadi
% Omidreza Amrollahinasab
% Montanuniversität Leoben, Austria
% Chair of Reservoir Engineering
% https://dpe.ac.at/
% ----------------------------------
%
%%
    figTitle   = data.title;
    yLabel     = data.yLabel;
    xLabel     = data.xLabel;
    figStyle   = data.style;
    figMarker  = data.marker;
    figRefresh = data.refresh;
    figTag     = data.tag;
    axisSide   = data.yyaxis;
    xLim       = data.xLim;
    yLim       = data.yLim;
    yScale     = data.yScale;
    if(figRefresh)
        ax = gca;
        hold(ax, "on");
    else
        % in case we want to create separate figures
%         fig = figure('Name',        figTitle, ...
%                      'Tag',         figTag,   ...
%                      'NumberTitle', 'off',     ...
%                      'WindowStyle', figStyle );
%         figure(fig);

        % for combined plots 
        nexttile
        ax = gca;
    end    
    if (axisSide == "right")
        yyaxis right
    end
    if strcmp(yScale,'linear')
        p = plot(ax,data.x,data.y,figMarker);
    elseif strcmp(yScale,'log')
        ax.YScale = 'log';
        p = semilogy(ax,data.x,data.y,figMarker);
    end
    ax.YColor = 'k';
    title(figTitle);
    xlabel(xLabel); ylabel(yLabel);
    if (~isempty(xLim))
        xlim(xLim);
    end
    if (~isempty(yLim))
        ylim(yLim);
    end
    grid on;
end