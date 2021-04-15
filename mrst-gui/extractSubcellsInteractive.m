function varargout = extractSubcellsInteractive(G, varargin)
%Extract a subset of cells from a grid interactively.
%
% SYNOPSIS:
%   cells = extractSubcellsInteractive(G)
%   cells = extractSubcellsInteractive(G, rock.poro)
%   extractSubcellsInteractive(G)
%
% DESCRIPTION:
%   Interactive extraction of a subset of cells.
%   First click defines i/j index to show subsets. The second plot can then
%   be clicked to toggle selection. Middle mouse can select entire column.
%
% REQUIRED PARAMETERS:
%   G - Valid grid structure.
%
% OPTIONAL PARAMETERS:
%  data - Plot a dataset onto the grid to make selection easier.
%
% RETURNS:
%   Selected subset.
%
% NOTE:
%   This function is under active development and is subject to change.
%
% EXAMPLE:
%   G = cartGrid([10,10,3])
%   G = computeGeometry(G);
%   c = extractSubcellsInteractive(G, rand(G.cells.num, 1))
%

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

    assert(isfield(G, 'cartDims'), 'Valid cartDims required');
    assert(all(isfinite(G.cartDims)), 'Valid cartDims required');
    assert(isfield(G.cells, 'centroids'), ['No geometry information found.'...
        ' Did you forget to call computeGeometry?']);
    hasOutput = nargout > 0;
    mainFig = figure;
    mainAx = NaN;
    reportFig = NaN;

    selectionFig = NaN;
    currentCell = NaN;

    if nargin > 1
        dataset = varargin{1};
        varargin = varargin(2 : end);
        plotG = @(G, varargin) plotCellData(G, dataset, varargin{:});
    else
        plotG = @(G, varargin) plotGrid(G, varargin{:});
    end


    opt = struct('selectedCells', [],...
                 'simple',        false, ...
                 'selcount',      inf);
    opt = merge_options(opt, varargin{:});

    bf = boundaryFaces(G);

    if ~isempty(opt.selectedCells)
        selectedCells = opt.selectedCells;
    else
        selectedCells = {[]};
    end


    ijk = cell(size(G.cartDims));
    [ijk{:}] = ind2sub(G.cartDims, G.cells.indexMap(1:G.cells.num));

    plotMain();

    if hasOutput
        uiwait(mainFig)
        if ishandle(reportFig), close(reportFig); end;
        if ishandle(selectionFig), close(selectionFig); end;
        varargout{1} = selectedCells;
    end

function plotMain()
    figure(mainFig); cla;
    mainAx = gca;
    plotG(G, 'ButtonDownFcn', @onClickPlot, 'FaceAlpha', .3);
    if ~isnan(currentCell)
        plotG(G, ijk{1} == ijk{1}(currentCell) | ijk{2} == ijk{2}(currentCell),...
                 'FaceColor', 'none', 'EdgeColor', 'White')
    end


    for i = 1:numel(selectedCells)
        subs = selectedCells{i};
        if any(subs)
            plotG(G, subs, 'EdgeColor', getColor(i), 'LineWidth', 2);
        end
    end
    axis tight off
    writeReport()
end

function writeReport()
    if isnan(currentCell) || hasOutput
        return
    end
    if ~ishandle(reportFig)
        reportFig = figure;
    end
    clf(reportFig);
    ph = uipanel(reportFig, 'Title', 'Selection info');

    s = sprintf('Current selection at %d, %d\n', (ijk{1}(currentCell)),...
                                                 (ijk{2}(currentCell)));

    for i = 1:numel(selectedCells)
        sc = selectedCells{i};
        if ~isempty(sc)
            s = [s 's' num2str(i) ' = [' sprintf('%d ', sc) ']' sprintf('\n')];
        end
    end

    uicontrol(ph, 'Units', 'Normalized', 'String', s, 'Style', 'edit', 'Position', [0 0 1 1], 'Max', 1e9)


end

function onClickPlot(src, event)
    pts = get(mainAx, 'CurrentPoint');
    [currentCell f] = nearestCellLine(G, bf, pts);

    if ~ishandle(selectionFig)
        selectionFig = figure;
    else
        figure(selectionFig); clf;
    end
    seli = ijk{1} == ijk{1}(currentCell);
    selj = ijk{2} == ijk{2}(currentCell);
    column = seli & selj;
    if opt.simple
        selectedCells = {find(column)};
    end

    subplot(2,1,1);
    plotSubset(seli, column)

    subplot(2,1,2);
    plotSubset(selj, column)

    plotMain();
end

function onClickSubset(src, event, subset)
    persistent currentSel;
    if isempty(currentSel); currentSel = 1; end;

    ax = get(src, 'Parent');
    f = get(ax, 'Parent');
    pts = get(ax, 'CurrentPoint');
    c = nearestCellLineSlice(G, subset, pts);

    seltype = get(f, 'SelectionType');

    if strcmpi(seltype, 'alt');
        if numel(selectedCells) < opt.selcount
            selectedCells = [selectedCells, []];
            currentSel = currentSel + 1;
        else
            currentSel = mod(currentSel, opt.selcount) + 1;
        end
    end
    if strcmpi(seltype, 'extend')
        c = find(ijk{1} == ijk{1}(c) & ijk{2} == ijk{2}(c));
    end
    alreadySelected = cellfun(@(x) any(ismember(c, x)), selectedCells);

    if ~any(alreadySelected)
        selectedCells{currentSel} = [reshape(selectedCells{currentSel}, 1, []),...
                                     reshape(c, 1, [])];
    else
        tmp = selectedCells{alreadySelected};
        tmp(ismember(tmp, c)) = [];
        selectedCells{alreadySelected} = tmp;
    end
    onClickPlot();
    plotMain();
end

function setView(selection)
    % Try to get a "nice" sideways view...
    bfs = boundaryFaces(G, selection);
    normals = sum(abs(G.faces.normals(bfs, :))./repmat(G.faces.areas(bfs), 1, 3), 1);
    % view uses angles, atan uses radians
    az = 360*atan(normals(1)/normals(2))/(2*pi);
    if isnan(az), az = 0; end;
    view(az, 0);
end

function plotSubset(selection, column)
    if islogical(selection)
        selection = find(selection);
    end

    handler = @(src, event) onClickSubset(src, event, selection);
    plotG(G, selection, 'ButtonDownFcn', handler, 'EdgeColor', 'k', 'EdgeAlpha', .25);
    plotG(G, column, 'EdgeColor', 'white', 'LineWidth', 2, 'FaceColor', 'None');
    for i = 1:numel(selectedCells)
        subs = selectedCells{i};
        subs = subs(ismember(subs, selection));
        if any(subs)
            plotG(G, subs, 'EdgeColor', getColor(i), 'LineWidth', 2, 'ButtonDownFcn', handler);
        end
    end

    setView(selection);
    axis tight off
end



end
function c = getColor(index)
    colors = {'y', 'm', 'c', 'r', 'g', 'b', 'w', 'k'};
    c = colors{mod(index, numel(colors)) + 1};
end

function c = nearestCellLineSlice(G, cells, pts)
    if islogical(cells)
        cells = find(cells);
    end
    gc = G.cells.centroids(cells,:);
    mgc = min(gc);
    Mgc = max(gc);

    N = numel(cells);
%     h = abs([diff(xlim) diff(ylim) diff(zlim)]);
%     h = Mgc - mgc;
    x1 = repmat(pts(1,:), N, 1);
    x2 = repmat(pts(2,:), N, 1);
    x0 = gc;
    % Find distance projected on line with a fudge factor for strange grids
    linedist = sum(cross(x0-x1, x0-x2 ,2).^2, 2);
    ptdist = sum((x0-x1).^2, 2);
    dist = sqrt(ptdist - linedist);

    [val ind] = min(dist); %#ok backwards compatability
    c = cells(ind);
end

