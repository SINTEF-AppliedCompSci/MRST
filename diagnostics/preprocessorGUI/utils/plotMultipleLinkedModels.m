function varargout = plotMultipleLinkedModels(G,vals,clims,patchCells,...
    modelNames,isLog,s3,wells,injColors,prodColors,plottitle,axview,cax,varargin)
%Undocumented Utility Function

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

opt = struct('camlight',          [], ...
             'outlineGrid',       [], ...
             'filterMinVal',      [], ...
             'filterMaxVal',      []);
[opt, ~] = merge_options(opt, varargin{:}); 

nM = numel(modelNames);


% Choose number of subplots
switch true
    case nM == 1
        frows = 1;
        fcols = 1;
    case nM ==2
        frows = 1;
        fcols = 2;
    case nM > 2 & nM <=4
        frows = 2;
        fcols = 2;
    case nM > 4 & nM <=6
        frows = 2;
        fcols = 3;
    case nM > 6 & nM <= 9   
        frows = 3;
        fcols = 3;
    case nM > 9 & nM <= 16
        frows = 4;
        fcols = 4;
    otherwise
        disp('Invalid number of models')
end

fh = figure('SizeChangedFcn',@resizeModelViewer);
fh.Name = plottitle;
axlist = [];

for i = 1:nM

    % Plot data
    d.Axes3D = subplot(frows,fcols,i);
    d.Patch = CellDataPatch(G, zeros(G.cells.num,1), ...
        'Parent', d.Axes3D, 'EdgeColor', [.3 .3 .3], ...
        'EdgeAlpha', .5, 'BackFaceLighting', 'lit');
  
    d.Patch.colorData = vals{i};
    d.Patch.cells = patchCells;  
    
    % Filter cells (code copied from filterPropCallback)
    if ~isempty(opt.filterMinVal)
        currentFilterRegion = and(opt.filterMinVal <= vals{i}, opt.filterMaxVal >= vals{i});
        d.Patch.cells = currentFilterRegion;
    end
    
    
    % Add outline grid
    if ~isempty(opt.outlineGrid)
        if strcmp(opt.outlineGrid.Visible,'on')
            plotGrid(G, 'FaceColor', 'none', 'EdgeAlpha', 0.15, 'EdgeColor', [.3 .3 .3]);
        end
    end
    
    d.WellPlot = WellPlotHandle(G, wells, ...
        'Visible', 'off', 'Parent', d.Axes3D);
    for k=1:numel(d.WellPlot.producers)
        d.WellPlot.producers(k).label.FontSize = 8;
    end
    for k=1:numel(d.WellPlot.injectors)
        d.WellPlot.injectors(k).label.FontSize = 8;
    end

    d.WellPlot.visibleInjectors = s3.wsel.injectorIx;
    d.WellPlot.visibleProducers = s3.wsel.producerIx;
    

    % Setup axes
    axis(d.Axes3D, 'tight', 'vis3d', 'off');
    d.Axes3D.ZDir = 'reverse';
    view(d.Axes3D, 3);
    daspect(d.Axes3D, [1 1 .2]);
    set(d.Axes3D,'View',axview);

    title(modelNames(i));

    % Add colorbar
    d.colorBar = colorbar(d.Axes3D,'WestOutside','AxisLocation', 'in');
    set(d.colorBar,'Units', 'normalized', 'Xaxis', 'left');
    d.colormap3D = d.Axes3D.Colormap;
    d.injColors  = injColors;
    d.prodColors = prodColors;
    
    if clims(1) ~= clims(2)
        updateColorBar(d, s3, clims, isLog);
    end
    set(d.Axes3D,'CLim',cax);
    % Add histogram
    if size(vals{i},2) == 1    
%         posHAx = d.colorBar.Position + [(sum(d.colorBar.Position([1 3])) + 0.002) 0 0 0];
        posHAx = d.colorBar.Position;

        d.colorHAx = axes('Position',posHAx, 'Units', 'normalized'); 
        d.colorHAx.Position(1) = sum( d.colorBar.Position([1 3])) + 0.002;        
        updateColorHist(d);
    end
    


    % Add lighting 
    if ~isempty(opt.camlight)
        if strcmp(opt.camlight.Visible,'on')
            camlight(d.Axes3D);
        end
    end
    

        
        
    
    ax{i} = d;
    axlist = [axlist,d.Axes3D];
    
    
    
end


Link = linkprop(axlist,{'CameraUpVector', 'CameraPosition', 'CameraTarget', ...
    'XLim', 'YLim', 'ZLim', 'CLim'});
setappdata(gcf, 'StoreTheLink', Link);


if nargout > 0
        varargout{1} = fh;
        if nargout > 1
            varargout{2} = ax;
        end
end


function updateColorHist(d)
    vals = d.Patch.colorData(d.Patch.cells);
    lims = d.Axes3D.CLim;
    [counts,bins] = hist(vals, linspace(lims(1), lims(2),25));
    visState = 'on';
    c = d.colorHAx.Children;
    if numel(c) == 1 && isprop(c, 'Visible')
        visState = d.colorHAx.Children.Visible;
    end
    barh(d.colorHAx, bins, counts,'hist');
%     d.colorHAx.XDir = 'reverse';
%     d.colorHAx.XAxisLocation = 'top';
    d.colorHAx.Tag = 'histogram';    
    axis(d.colorHAx,'tight', 'off');
    h = findobj(d.colorHAx,'Type','patch');
    set(h,'FaceColor','none','EdgeColor',[.2 .2 .2], 'Visible', visState);

    
function updateColorBar(d, s3, lims, islog)

    cb = d.colorBar;
    
    % let first matlab choose ticks/labels
    if ~(s3.psel.typeIx == 3 && any(s3.psel.propIx == [7,8]))
        cb.TicksMode = 'auto';
        cb.TickLabelsMode = 'auto';
        cb.Limits = lims;
        if islog
            p = cb.Ticks;
            cb.TickLabels = cellfun(@(x)sprintf('%2.1e',10^x), num2cell(p), ...
                'UniformOutput', false);
            cb.TicksMode = 'manual';
        end
%         set(d.colorHAx.Children, 'Visible', 'on');
%         d.setColormap3D('default');
        %cb.Position(1) = d.layoutParams.menuWidth + 50;
        %d.colorHAx.Position(1) = d.colorBar.Position(1)+d.colorBar.Position(3)+10;
    else
        if  s3.psel.propIx == 7 % sweep regions
            ninj = numel(d.WellPlot.injectors);
            setColormap3D(d,'injectors');
            cb.Ticks = (.5:ninj)/ninj;
            cb.TickLabels = s3.wsel.injSelector.String;
            cb.Limits = [0 1];
        elseif s3.psel.propIx == 8 % drainage regions
            nprod = numel(d.WellPlot.producers);
            setColormap3D(d,'producers');
            cb.Ticks = (.5:nprod)/nprod;
            cb.TickLabels = s3.wsel.prodSelector.String;
            cb.Limits = [0 1];
        end
%         set(d.colorHAx.Children, 'Visible', 'off');
        cb.TickLabelInterpreter = 'none';
        %cb.Position(1) = d.layoutParams.menuWidth + 75;
        % switch off log if selected
        if s3.psel.logSwitch == true
            s3.psel.logSwitch = false;
            s3.psel.logSwitchCallback(s3.psel.logSwitchBox, []);
        end
    end
    d.colorBar = cb;
    
 function setColormap3D(d, str)
        if ~strcmp(str, d.colormap3D)
            if strcmp(str, 'injectors')
                ninj = numel(d.WellPlot.injectors);
                cmap = d.injColors(1:ninj,:);
                colormap(d.Axes3D, cmap);
            elseif strcmp(str, 'producers')
                nprod = numel(d.WellPlot.producers);
                cmap = d.prodColors(1:nprod,:);
                colormap(d.Axes3D, cmap);
            else
                colormap(d.Axes3D, str);
            end
            d.colormap3D = str;
        end



function [lims, vals, flag] = makeSafeForLog(lims, vals, oom)
if lims(2) <= 0
    flag = false;
else
    flag = true;
    if lims(1) <= 0 || lims(1) < lims(2)*10^(-oom)
        lims(1) = lims(2)*10^(-oom);
        vals = max(vals, lims(1));
    end
end
