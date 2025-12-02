function PlotKrPc(model)
%
% DESCRIPTION: plot the relative permeability and capillary pressure
%
% SYNOPSIS:
%   PlotKrPc(model)
%
% PARAMETERS:
%   model - struct containing following fields:
%   - experiment: saturation functions used for forward modeling
%
% RETURNS:
%   plot of the capillary pressure and relative permeability
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
satfun = model.satfun;

% create a tiled figure
fig = figure('Name', 'Saturation Functions used for forward simulation', ...
            'NumberTitle', 'off', ...
            'WindowStyle', 'docked' );
tiledlayout(2,2)

%% plot relative permeabilities
x  = satfun.sw_kr;
y1  = satfun.kro;
y2  = satfun.krw;
headers = {'Sw [fraction]','Kro [fraction]','Krw [fraction]'};    
[headers,units] = SplitHeaders(headers);
data.xLabel = strcat(headers{:,1},32,units{:,1});
data.yLabel = strcat(headers{:,2},32,units{:,2});
yyLabel      = strcat(headers{:,3},32,units{:,2});    
data.title  = 'Relative Permeability (Linear scale)';
data.style = 'docked';
data.refresh = 0;
data.tag     = data.title;

% linear scale
data.yScale = 'linear';
data.xLim    = [0 1];
data.yLim    = [0 1];
data.yyaxis  = "left";
data.x = x;
data.y = y1;
data.marker  = '-r.';
PlotData(data);
data.refresh = 1;
data.yyaxis  = "right";
data.yLabel  = yyLabel;
data.y = y2;
data.marker  = '-b.';
PlotData(data);

% log scale
data.title  = 'Relative Permeability (Log Scale)';
data.refresh = 0;
data.tag     = data.title;
data.yScale = 'log';
data.xLim    = [0 1];
data.yLim    = [1e-5 1];
data.yLabel = strcat(headers{:,2},32,units{:,2});
data.yyaxis  = "left";
data.y = y1;
data.marker  = '-r.';
PlotData(data);
data.refresh = 1;
data.yyaxis  = "right";
data.yLabel  = yyLabel;
data.y = y2;
data.marker  = '-b.';
PlotData(data);

%% plot capillary pressure

% linear plot
data.x = satfun.sw_pc;
data.y = satfun.pc/1e5;    
headers = {'Sw [fraction]','Pc [bar]'};    
[headers,units] = SplitHeaders(headers);
data.xLabel  = strcat(headers{:,1},32,units{:,1});
data.yLabel  = strcat(headers{:,2},32,units{:,2});   
data.title   = 'Capillary Pressure';
data.style = 'docked';
data.marker  = '-b.';
data.refresh = 0;
data.tag     = data.title;
data.yyaxis  = "left";
data.xLim    = [0 1];
if min(data.y) >= 0
    data.yLim    = [min(data.y)-0.1 min(max(data.y),2) + 0.1];
else
    data.yLim    = [max(min(data.y),-2) min(max(data.y),2)];
end
data.yScale = 'linear';
PlotData(data);    

% log plot
data.yScale = 'log';
data.title   = 'Capillary Pressure (Log Scale)';
if min(data.y) >= 0
    data.yLim    = [0 min(max(data.y),2)+0.1];
else
    data.yLim    = [0 min(max(data.y),2)];
end
if not(any(data.y<0))
    PlotData(data); 
end