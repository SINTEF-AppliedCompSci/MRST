function PlotData(data)
% <keywords>
%
% Purpose : plot data
%
% Syntax :
%   PlotData(data)
%
% Input Parameters :
%   data: struct containing arrays to be plotted and the plotting features
%
% Return Parameters :
%   
%   
% Description :
%
% Author : 
%    Siroos Azizmohammadi
%    Omidreza Amrollahinasab
%
% History :
% \change{1.0}{09-Nov-2021}{Original}
%
% --------------------------------------------------
% (c) 2021, Siroos Azizmohammadi,
% Omidreza Amrollahinasab
% Chair of Reservoir Engineering, University of Leoben, Austria
% email: info@dpe.ac.at
% url: dpe.ac.at
% --------------------------------------------------
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