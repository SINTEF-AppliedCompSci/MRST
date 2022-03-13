function h = mapPlot(h, Gt, varargin)
% Plot a map of the 2D formation, with optional features.
% NB: Make sure to open a figure and set 'hold on' before calling this function.
%
% SYNOPSIS:
%   function mapPlot(Gt, varargin)
%
% DESCRIPTION:
%
% PARAMETERS:
%   Gt       - 
%   varargin - 
%
% RETURNS:
%   mapPlot(Gt, - 

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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

    opt.traps = [];
    opt.trapsets = []; % if we want multiple colors on traps, this should be
                       % a cell array with trap indices for each color
                       % group.  (Default is a single set).
    opt.trapcolor = [1 0 0]; % if we want multiple colors on traps, there
                             % should be one line per trapset. 
    opt.trapalpha = 0.3;
    opt.rivers = []; % 'cell_lines' in trap structure
    opt.rivercolor = [1 0 0];
    opt.plumes = [];
    opt.labels = [];
    opt.title = [];
    opt.wellcells = [];
    opt.well_numbering = false;
    opt.plume_h_threshold = 0.3;
    opt.plumeLineWidth = 2;         % original value passed into drawContour for 'lineWidth'
    opt.plumeLineStyle = '-';       % default of contour() used by drawContour()
    opt.casecolors = 'rbcgm';
    opt.maplines = 40;
    opt.quickclip = true;
    opt.background = [];
    opt.backgroundalpha = 1;
    opt.background_threshold = [];
    opt.colorbar = 'background';
    opt.colorbarposition = 'South';
    
    opt = merge_options(opt, varargin{:});
    figure(h); hold on;
    ax = get(h, 'currentaxes');
    
    % Setting title
    if(~isempty(opt.title))
        title(opt.title,'FontSize', 14);
    end
    % plotting map
    drawContours(ax, Gt, Gt.cells.z, opt.maplines, 'color', 'k');
    
    % drawing background
    if ~isempty(opt.background)
       % eliminating any negative infinite values
       neginf = opt.background == -Inf;
       opt.background(neginf) = min(opt.background(~neginf));
       ax = drawSmoothField(ax, Gt, ...
                            opt.background , 200           , ...
                            'quickclip'    , opt.quickclip , ...
                            'alpha'        , opt.backgroundalpha, ...
                            'field_threshold', opt.background_threshold);
       if strcmpi(opt.colorbar, 'background')
           caxis([min(opt.background), max(opt.background)]);
           hh = colorbar('peer', ax, opt.colorbarposition);
           set(hh, 'fontSize', 14);
       end
    end    
    
    % Plotting plumes, if any
    num_plumes = size(opt.plumes, 2);
    if num_plumes > numel(opt.casecolors)
        error('Not enough colors to draw %i plumes', num_plumes);
    end

    for c_ix = 1:num_plumes
        drawContours(ax, ...
                     Gt, opt.plumes(:, c_ix), ...
                     [opt.plume_h_threshold, opt.plume_h_threshold], ...
                     'color', opt.casecolors(c_ix), ...
                     'lineWidth', opt.plumeLineWidth, ...
                     'lineStyle', opt.plumeLineStyle, ...
                     'quickclip', opt.quickclip);
    end
    
    % Plotting plume labels
    if num_plumes > 1
        legend(opt.labels, 'Location', 'SouthEast');    
    end
    
    % plotting wellcells
    plot(Gt.cells.centroids(opt.wellcells,1), ...
         Gt.cells.centroids(opt.wellcells,2), ...
         'ok', 'MarkerSize',10, 'MarkerFaceColor', 'k','MarkerEdgeColor', ...
         'y');
    if opt.well_numbering
        xpos = Gt.cells.centroids(opt.wellcells,1);
        ypos = Gt.cells.centroids(opt.wellcells,2);
        labels = [repmat(' ', numel(xpos), 1), num2str([1:numel(xpos)]')]; %#ok
        text(xpos, ypos, labels, 'fontsize', 16);
    end
    
    % Plotting traps
    if ~isempty(opt.traps)
        % colorize traps
        if isempty(opt.trapsets) || size(opt.trapcolor, 1) == 1
            % all traps will have same color
            colorizeRegion(h, Gt, opt.traps, opt.trapcolor(1,:), 'faceAlpha' , opt.trapalpha);
        else
            % check that we have enough colors
            assert(size(opt.trapcolor,1) >= numel(opt.trapsets)); 
            color_ix = 1;
            for ts = opt.trapsets
                traps = ismember(opt.traps, [ts{:}]);
                colorizeRegion(h, Gt, traps, opt.trapcolor(color_ix,:), ...
                               'faceAlpha', opt.trapalpha);
                color_ix = color_ix + 1;
            end
        end
    end
    if ~isempty(opt.rivers)
        % sketch rivers
        for tr = opt.rivers
            for r = tr{:}
                draw_cell_connections(Gt, r{:}, 'color', opt.rivercolor, 'lineWidth', 1);
            end
        end
    end
    

    
end

% ============================================================================

% ----------------------------------------------------------------------------
function draw_cell_connections(Gt, cells, varargin)

    x = Gt.cells.centroids(cells, 1);
    y = Gt.cells.centroids(cells, 2);
    smooth=@(x) x;
    x = smooth(x);
    y = smooth(y);
    plot(x, y, varargin{:});
    
end
