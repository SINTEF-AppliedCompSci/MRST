 function varargout = plotToolbar(Grid, inputdata, varargin)
%Plot one or more datasets with interactive tools for visualization
%
% SYNOPSIS:
%   h = plotToolbar(G, dataset)
%   h = plotToolbar(G, struct('data1', data, 'data2', data));
%
% DESCRIPTION:
%   plotToolbar is designed to be a drop-in replacement for plotCellData
%   which has several added features that allow the user to inspect and
%   manipulate several datasets. The key features are:
%
%       - Capable of interacting with multiple datasets simultanously. See
%       required parameters.
%
%       - Supports all options supported by plotCellData.
%
%       - Improved rotation speed of figures due to custom rotate handler.
%       This is extremely useful when working with 100k+ cells.
%
%       - Subdivide model, use slicing planes to visualize internal
%       features or limit by logical indices or the values of datasets.
%
%       - Dynamic histograms allow easy inspection of subsamples of a
%       complex dataset.
%
%       - Multiple transformations (log10, absolute value etc) makes it
%       easy to explore a new dataset without any a priori knowledge of the
%       data. Visualize directly on the equivialent logical grid structure.
%
%       - Visualize a whole simulation, including playback with dynamic
%       selection.
%
% REQUIRED PARAMETERS:
%   Grid    - Grid data structure.
%
%   data    - Scalar cell data with which to colour the grid.  One scalar,
%             indexed colour value for each cell in the grid or one
%             TrueColor value for each cell.  If a cell subset is specified
%             in terms of the 'cells' parameter, 'data' must either contain
%             one scalar value for each cell in the model or one scalar
%             value for each cell in this subset.
%
%   cells   - Vector of cell indices defining sub grid.  The graphical
%             output of function 'plotToolbar' will be restricted to the
%             subset of cells from 'G' represented by 'cells'.
%
%             If unspecified, function 'plotToolbar' will behave as if the
%             user defined
%
%                 cells = 1 : G.cells.num
%
%             meaning graphical output will be produced for all cells in
%             the grid model 'G'.  If 'cells' is empty (i.e., if
%             ISEMPTY(cells)), then no graphical output will be produced.
%
%   'pn'/pv - List of property names/property values.  OPTIONAL.
%             This list will be passed directly on to function PATCH
%             meaning all properties supported by PATCH are valid.
%
% OPTIONAL PARAMETERS:
%   
%   Most supported options correspond to activating buttons in the GUI
%   at the first launch. This is useful when the user wants to specify the
%   exact setup for the plotting window.
%
%   log10  - Toggle default logarithm transform of data.
%
%   exp    - Toggle default exponential transform of data.
%
%   abs    - Toggle default absolute value transform of data.
%
%   filterzero - Toggle filtering of near-zero values.
%
%   outline - Toggle default outline of grid.
%
%   pauseTime - Time to use between plot updates when interactive playback
%               is used.
%
%   lockCaxis - Toggle locking of caxis for dynamic playback.
%
%   plot1d    - Toggle 1D plotting (requires grid with only one dimension
%               with more than one cell thickness)
%
%   plotMarkers - Markers used for 1D plotting. Defaults to none.
%
%   field       - Initial field to plot. For multiple columns, the syntax 
%                 's:3' is supported (take column 3 of matrix s). Must be a
%                 valid plotting field.
%
%   startplayback - Start playback immediately with no user input.
%
% RETURNS:
%   h       - Handle to underlying plot call. If used interactively, this
%             handle will *not* remain usable as this static value will not
%             be updated when the plot changes.
%
% EXAMPLE:
%    G = cartGrid([10,10,10])
%    G = computeGeometry(G);
%    %Visualize the grid itself
%    plotToolbar(G, G)
%
% SEE ALSO:
%   `plotCellData`, `plotFaces`

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
    G = Grid;
    if mod(numel(varargin), 2) == 1
        initialSelection = varargin{1};
        if ~islogical(initialSelection)
            tmpsel = false(G.cells.num, 1);
            tmpsel(initialSelection) = true;
            initialSelection = tmpsel;
            clear tmpsel
        end
        assert(numel(initialSelection) == G.cells.num);
        initialSelection = reshape(initialSelection, [], 1);
        varargin = varargin(2 : end);
    else
        initialSelection = true(G.cells.num, 1);
    end
    opt = struct('log10',         false, ...
                 'exp',           false, ...
                 'abs',           false, ...
                 'filterzero',    false, ...
                 'time',          [],    ...
                 'dynamicTitle',  false, ...
                 'title',         '',    ...
                 'logical',       false, ...
                 'outline',       false, ...
                 'pauseTime',     0.150, ...
                 'lockCaxis',     false, ...
                 'plot1d',        false, ...
                 'plotmarkers',   false, ...
                 'skipAugmented', false, ...
                 'field',         '',    ...
                 'step_index',    1,     ...
                 'startplayback', false  ...
                 );
             
    [opt, varargin] = merge_options(opt, varargin{:});

    isscalar = false;
    
    N = numel(inputdata);
    accessdata = @(x) inputdata(x);
    if isPlottable(inputdata) && max(size(inputdata)) == G.cells.num
        if size(inputdata, 2) == 1
            isscalar = true;
        end
        inputdata = struct('values', inputdata);
        N = 1;
        accessdata = @(x) inputdata;
    elseif iscell(inputdata)
        if all(cellfun(@isPlottable, inputdata))
            % Transform into structs
            inputdata = cellfun(@(x) struct('data', x), inputdata, 'UniformOutput', false);
        end
        accessdata = @(x) inputdata{x};
    elseif strcmpi(class(inputdata), 'ResultHandler')
        N = inputdata.numelData();
        accessdata = @(x) inputdata{x};
    end

    datasetname = inputname(2);
    [cellfields, hasCell, hasNode] = getStructFields(G, accessdata(1), datasetname); %#ok
    if (hasNode && isfield(G, 'nodes')) && ~isfield(G.nodes, 'cells') && ~opt.skipAugmented
        try
            G = createAugmentedGrid(G);
        catch
            disp('Unable to create augmented grid, node filtering not available');
        end
    end
    ni = min(opt.step_index, N);
    if isempty(cellfields)
        error('Input dataset contains no fields with G.cells.num rows suitable for plotting.');
    else
        if ~isempty(opt.field)
            ix = find(strcmpi(cellfields, [datasetname, '.', opt.field]));
            if isempty(ix)
                ix = 1;
            end
        else
            ix = 1;
        end
        selectionName = cellfields{ix};
    end
    data = readStructField(accessdata(ni), selectionName);
    % Constants
    if ~isfield(G, 'cartDims')
        G.cartDims = ones(1, G.griddim);
    end
    
    if G.griddim == 2
        if isfield(G.faces, 'nodes')
            G.cells.sortedCellNodes = getSortedCellNodes(G);
        elseif isfield(G, 'parent') && isfield(G.parent, 'nodes')
            G.parent.cells.sortedCellNodes = getSortedCellNodes(G.parent);
        end
    end
    
    ijk = cell(size(G.cartDims));
    if isfield(G.cells, 'indexMap')
        [ijk{:}] = ind2sub(G.cartDims, G.cells.indexMap(1:G.cells.num));
    end
    % Plotting parameters
    [...
     vecToggle, ...
     linePlotToggle, ...
     markerToggle, ...
     logDisplay, ...
     absDisplay, ...
     tenDisplay, ...
     nonzeroDisplay, ...
     slices, ...
     logicalsubset, ...
     histsubset, ...
     histh, ...
     cax, ...
     histAll, ...
     N_hist, ...
     sliceoutlineh, ...
     gridOutline, ...
     minmax, ...
     mainFig, ...
     ph, ...
     vh, ...
     G_logi...
     ] = deal([]);
    
    fig = gcf;
    % Delete previous instances of toolbar in same figure
    delete(findobj(fig, 'Tag', 'mrst-plottoolbar'));

    ht = uitoolbar(fig, 'Tag', 'mrst-plottoolbar');

    % Log10 button
    log10toggleh = uitoggletool(ht, 'TooltipString', 'Logarithmic plot (log10)',...
                     'ClickedCallback', @toggleTransform, ...
                     'cdata', geticon('log10'), 'State', boolToOnOff(opt.log10));
    % 10^x button
    tentoggleh = uitoggletool(ht, 'TooltipString', '10^x plot',...
                     'ClickedCallback', @toggleTransform, ...
                     'cdata', geticon('10'), 'State', boolToOnOff(opt.exp));

    % Abs button
    abstoggleh = uitoggletool(ht, 'TooltipString', 'Absolute value',...
                     'ClickedCallback', @toggleTransform, ...
                     'cdata', geticon('abs'), 'State', boolToOnOff(opt.abs));

    % Abs button
    nonzerotoggleh = uitoggletool(ht, 'TooltipString', 'Filter zero values',...
                     'ClickedCallback', @toggleTransform, ...
                     'cdata', geticon('nonzero'), 'State', boolToOnOff(opt.filterzero));

    % Logical grid button
    uitoggletool(ht, 'TooltipString', 'Logical grid',...
                     'ClickedCallback', @toggleLogicalGrid, ...
                     'cdata', geticon('ijkgrid'), 'State', boolToOnOff(opt.logical));

    % Grid outline button
    gridToggle = uitoggletool(ht, 'TooltipString', 'Show grid boundary outline',...
                   'ClickedCallback', @replotPatch,...
                   'cdata', geticon('outline'), 'State', boolToOnOff(opt.outline));

    % Histogram button
    uipushtool(ht, 'TooltipString', 'Dynamic histogram of currently displayed data',...
                   'ClickedCallback', @plotHistogram,...
                   'cdata', geticon('hist'));
    % Limiting dataset
    uipushtool(ht, 'TooltipString', 'Set limits on displayed dataset',...
                   'ClickedCallback', @limitDataset,...
                   'cdata', geticon('minmax'),...
                   'Separator', 'on');
    % Limit logical
    uipushtool(ht, 'TooltipString', 'Pick logical indices',...
                   'ClickedCallback', @limitLogical,...
                   'cdata', geticon('ijk'));

    % Slice!
    sliceToggle = uitoggletool(ht, 'TooltipString', 'Slice the figure. Click two times to slice.',...
                   'ClickedCallback', @sliceToggleFn,...
                   'cdata', geticon('slice'));

    % Slice 2!
    uipushtool(ht, 'TooltipString', 'Interactive plane slice',...
                   'ClickedCallback', @(src, event) plotAdjustiblePlane(G, data, 'Callback', @planarCallback),...
                   'cdata', geticon('sliceplane'));

    % Reset selection
    uipushtool(ht, 'TooltipString', 'Reset selection',...
                   'ClickedCallback', @(src, event) resetSelection,...
                   'cdata', geticon('reset'));

    % Adjust patch
    uipushtool(ht, 'TooltipString', 'Set patch properties',...
                   'ClickedCallback', @(src, event) adjustPatch,...
                   'cdata', geticon('adjust'), ...
                   'Separator', 'on');


    % Axis adjustments
    uipushtool(ht, 'TooltipString', 'Set axis to tight',...
                   'ClickedCallback', @(src, event) axis('tight'),...
                   'cdata', geticon('tight'), ...
                   'Separator', 'off');

    % Zoom to selection
    zoomselToggle = uitoggletool(ht, 'TooltipString', 'Center axis on selection',...
                   'ClickedCallback', @(src, event) replotPatch(),...
                   'cdata', geticon('centersel'), ...
                   'Separator', 'off');

    % Freeze caxis
    caxisToggle = uitoggletool(ht, 'TooltipString', 'Freeze caxis',...
                   'cdata', geticon('lockcaxis'), ...
                   'Separator', 'off', 'State', boolToOnOff(opt.lockCaxis));

    if ~isfield(G, 'cartDims') || nnz(G.cartDims > 1) < 2 || opt.plot1d
        linePlotToggle = uitoggletool(ht, 'TooltipString', 'View as line plot',...
                       'ClickedCallback', @(src, event) replotPatch(),...
                       'cdata', geticon('1d'), ...
                       'Separator', 'off', 'State', boolToOnOff(opt.plot1d));
        markerToggle = uitoggletool(ht, 'TooltipString', 'Show markers in 1D plot',...
                       'ClickedCallback', @(src, event) replotPatch(),...
                       'cdata', geticon('marker'), ...
                       'Separator', 'off', 'State', boolToOnOff(opt.plotmarkers));
    end
    addSelector(ix);

    % Initial plot comes here
    resetSelection(varargin{:})

    if isnumeric(fig)
        % Override default drawing button to make plot rotations much faster
        % This is only needed for HG1.
        fastRotateButton();
    end

    if nargout > 0
        varargout{1} = ph;
        if nargout > 1
            varargout{2} = @(newGrid, newDataset, varargin) injectDataset(newGrid, newDataset, varargin{:});
        end
    end

    function closeFcn(src, event) %#ok<*INUSD>
        % Destructor
        delete(fig);
    end
    set(fig, 'CloseRequestFcn', @closeFcn);
    if exist('disableDefaultInteractivity', 'file')
        % Disable tooltips since we are doing 3D plotting without tooltips
        disableDefaultInteractivity(gca)
    end

    if opt.startplayback
        playbutton = findobj(fig, 'Tag', 'mrst-datasetselector-play');
        if ~isempty(playbutton)
            clickh = get(playbutton, 'Callback');
            clickh([], []);
        end
    end
 function injectDataset(newGrid, newDataset, varargin)
     G = newGrid;
     if numel(varargin)
         ni = varargin{1};
     end
     inputdata = newDataset;
     if iscell(inputdata)
         d = inputdata{ni};
     else
         d = inputdata(ni);
     end
     data = readStructField(d, selectionName);
     
     oldsel = findobj(fig, 'Tag', 'mrst-activedataset');
     addSelector(get(oldsel, 'value'));
     replotPatch();
     drawnow
 end

function planarCallback(selection, outlinepts)
    set(gridToggle, 'State', 'on')
    slices = selection;
    deleteHandle(sliceoutlineh);
    set(0, 'CurrentFigure', mainFig);
    axis tight off
    sliceoutlineh = patch(outlinepts(:, 1), outlinepts(:, 2), outlinepts(:, 3), 'k', 'FaceAlpha', 0);
    replotPatch();
end

function selectdata(datas, index, fullname)
    data = datas;
    if ~strcmpi(selectionName, fullname)
        set(caxisToggle, 'State', 'off');
        selectionName = fullname;
    end
    ni = index;

    replotPatch();
end

function h = getHandle()
    h = ph;
end

function addSelector(ix)
    if nargin == 0
        ix = 1;
    end
    if N > 1 || ~isscalar
        set(0, 'CurrentFigure', fig)

        oldPanel = findobj(fig, 'Tag', 'mrst-datasetselector');
        if ~isempty(oldPanel)
            delete(oldPanel)
        end
        dsheight = .06;
        hmod = .1;
        subplot('Position', [hmod, dsheight + hmod, 1 - 2*hmod, 1-dsheight - 2*hmod]);
        datasetSelector(G, inputdata, ...
            'Parent', gcf, ...
            'Location', [0, 0, 1, dsheight], ...
            'Callback', @selectdata, ...
            'Setname', datasetname, ...
            'active', ix, ...
            'N', ni, ...
            'pauseTime', opt.pauseTime);
        set(fig, 'ResizeFcn', @handleFigureResize);
        handleFigureResize(fig, []);
   end
end

function sliceToggleFn(src, event)
     activateuimode(gcf, '');
     switch get(src, 'State')
         case 'on'
            set(getHandle(), 'ButtonDownFcn', @sliceDataset)
            setptr(fig, 'add')
         case 'off'
            set(getHandle(), 'ButtonDownFcn', [])
            setptr(fig,'arrow');
     end
 end

function toggleLogicalGrid(src, event)
    if strcmpi(get(src, 'State'), 'on')
        if ~isstruct(G_logi)
            G_logi = reorderLogical(G);
        end
        G = G_logi;
    else
        G = Grid;
    end
    replotPatch();
    axis tight
end

function G2 = reorderLogical(G)
    G2 = cartGrid(G.cartDims);
    cellInd = false(G.cells.num, 1);
    cellInd(G.cells.indexMap) = true;
    G2 = extractSubgrid(G2, find(cellInd));
    try
        G2 = mcomputeGeometry(G2);
    catch
        G2 = computeGeometry(G2);
    end
    G2.cartDims = G.cartDims;
end

function resetSelection(varargin)
    logDisplay = strcmpi(get(log10toggleh, 'State'), 'on');
    absDisplay = strcmpi(get(abstoggleh, 'State'), 'on');
    tenDisplay = strcmpi(get(tentoggleh, 'State'), 'on');
    nonzeroDisplay = strcmpi(get(nonzerotoggleh, 'State'), 'on');

    slices = NaN;
    logicalsubset.i = NaN;
    logicalsubset.j = NaN;
    logicalsubset.k = NaN;
    minmax.filter = [];
    minmax.active = true(G.cells.num, 1);
    mainFig = gcf;
    histsubset = nan;
    
    deleteHandle(sliceoutlineh);
    
    replotPatch(varargin{:});
end

function toggleTransform(src, event)
    logDisplay = strcmpi(get(log10toggleh, 'State'), 'on');
    absDisplay = strcmpi(get(abstoggleh, 'State'), 'on');
    tenDisplay = strcmpi(get(tentoggleh, 'State'), 'on');
    nonzeroDisplay = strcmpi(get(nonzerotoggleh, 'State'), 'on');

    set(caxisToggle, 'State', 'off')
    replotPatch();
end

function plotHistogram(varargin)
    if isempty(histAll); histAll = true; end
    if isempty(N_hist); N_hist = 10; end
    oldfig = gcf;

    if isempty(histh) || ~ishandle(histh)
        if nargin == 0
            return;
        end
        histh = figure;
    end
    [d, fun, subset] = getData(); %#ok<ASGLU>
    d = double(d);
    if histAll
        data_hist = d;
    else
        data_hist = d(subset, :);
    end
    data_hist = data_hist(all(isfinite(data_hist), 2), :);
    % If vector data, assume l2 norm
    if size(data_hist, 2) > 1
        data_hist = sum(data_hist, 2).^2;
    end

    set(0, 'CurrentFigure', histh);
    clf;
    % Find data
    [n, X] = hist(data_hist, N_hist); %#ok

    function h = plotRectangle(x, height, width, value, varargin)
        h = patch([x; x+width; x+width; x], [0; 0; height; height], value, varargin{:});
    end

    h = X(2) - X(1);
    for i = 1:N_hist
        plotRectangle(X(i) - h, n(i), h, X(i), 'ButtonDownFcn', @(src, event) histClickHandler(src, event, X(i), h, d));
    end

    axis tight
    caxis(cax);


    function histToggleAll(src, event)
        histAll = get(src, 'Value');
        plotHistogram();
    end

    function histChangeN(src, event)
        strings = get(src, 'String');
        N_hist = str2double(strings{get(src, 'Value')});
        plotHistogram();
    end
    bins = [10, 25, 50, 100, 250, 1000];
    uicontrol(histh, 'Style', 'popupmenu', 'String', arrayfun(@(x) num2str(x), bins , 'UniformOutput', false),...
        'Position', [5 0 60 50], 'Callback', @histChangeN, 'value', find(bins == N_hist))
    uicontrol(histh, 'Style', 'togglebutton', 'String', 'All',...
        'Position', [5 55 60 25], 'Value', histAll, 'Callback', @histToggleAll)
    set(0, 'CurrentFigure', oldfig)

    function closeFcn(src, event)
        histh = NaN;
        delete(src);
    end

    set(histh, 'CloseRequestFcn', @closeFcn);
end


function histClickHandler(src, event, x, h, d) %#ok<INUSL>
    seltype = get(histh, 'SelectionType');
    inside = abs(d - x) <= 1.0001*h/2;
    switch seltype
        case 'normal'
            % Normal click
            histsubset = inside;
        case 'alt'
            % Right click / ctrl-click
            if isnan(histsubset)
                histsubset = inside;
            else
                histsubset = histsubset | inside;
            end
        case 'extend'
            % Shift click
            histsubset = ~inside;
    end
    replotPatch();
 end

function sliceDataset(src, event)
    % Two mouse events.
    % The first defines a line through the axis. The second the last point
    % for a slicing plane.
    % If the first click is alternate (ctrl-click / mouse2) it is additive.
    % If the second click is extend the slice does not terminate

    axold = axis();

    pts1 = get(gca, 'CurrentPoint');
    additive = ~strcmpi(get(gcf, 'SelectionType'), 'normal') & ~isnan(slices);
    hold on
    while 1
        while 1
            if waitforbuttonpress == 1
                return
            end
            final = ~strcmpi(get(gcf, 'SelectionType'), 'extend');
            pts2 = get(gca, 'CurrentPoint');
            break
        end

        A = pts1(1,:);
        B = pts2(1,:);
        C = pts1(2,:);

        Normal = cross(B-A, C-A);

        if additive
            old = slices;
        else
            old = true(G.cells.num, 1);
        end
        slices = old & sum(repmat(Normal, G.cells.num, 1).*(G.cells.centroids  - repmat(A, G.cells.num,1)), 2) > 0;
        replotPatch();
        axis(axold);
        if final
            setptr(fig,'arrow');
            set(sliceToggle, 'State', 'off')
            break;
        end
    end
end

function limitLogical(src, event)

    pos = get(gcf, 'OuterPosition');
    fi = figure('Position',[pos(1:2)-[300 0],  [300 340]], 'Toolbar','none', 'MenuBar', 'none');
    set(fi, 'Name', 'Select logical indices');


    Mv = max([ijk{:}]);

    f = {'i', 'j', 'k'};
    for i = 1:G.griddim
        if isnan(logicalsubset.(f{i}))
            logicalsubset.(f{i}) = 1:Mv(i);
        end
    end

    integer_labels = @(arr) arrayfun(@(i) sprintf('%d', i), arr, 'UniformOutput', false);

    si = uicontrol(fi, 'Position', [0,   40, 100, 300], 'String', integer_labels(1:Mv(1)) , 'Style', 'listbox', 'Max', Mv(1), 'Value', logicalsubset.i, 'TooltipString', 'I');
    sj = uicontrol(fi, 'Position', [100, 40, 100, 300], 'String', integer_labels(1:Mv(2)) , 'Style', 'listbox', 'Max', Mv(2), 'Value', logicalsubset.j, 'TooltipString', 'J');
    if G.griddim == 3
        sk = uicontrol(fi, 'Position', [200, 40, 100, 300], 'String', integer_labels(1:Mv(3)) , 'Style', 'listbox', 'Max', Mv(3), 'Value', logicalsubset.k, 'TooltipString', 'K');
    end

    uicontrol(fi, 'Style', 'pushbutton', 'Position', [100, 0, 100, 40], 'string', 'Apply', 'callback', @(src, event) uiresume(fi))
    uiwait(fi);

    if ~ishandle(fi)
        % User closed window
        return
    end

    logicalsubset.i = get(si, 'Value');
    logicalsubset.j = get(sj, 'Value');
    if G.griddim == 3
    logicalsubset.k = get(sk, 'Value');
    end

    replotPatch();
    close(fi);
    limitLogical(src, event);
end

function limitDataset(src, event)
    addFilters(G, minmax.filter, accessdata(ni), @setlimits);
end

function setlimits(active, filter)
    % seperate function to break closure
    minmax.filter = filter;
    minmax.active = active;
    replotPatch();
end

function adjustPatch(src, event)
    % New figure
    pos = get(gcf, 'OuterPosition');
    fi = figure('Position',[pos(1:2)-50,  [350 120]], 'Toolbar','none', 'MenuBar', 'none');

    set(gcf, 'Name', 'Change patch style');

    se = linkedSlider(fi, [10,60], 0, 1, get(ph, 'EdgeAlpha'), 'Edge');
    sf = linkedSlider(fi, [10,90], 0, 1, get(ph, 'FaceAlpha'), 'Face');
    applyfun = @(src, event) set(ph, 'FaceAlpha', get(sf, 'Value'), 'EdgeAlpha', get(se, 'Value'));

    patchcolor = @(src, event) set(ph, 'EdgeColor', uisetcolor(get(ph, 'EdgeColor')));
    uicontrol(fi, 'Style', 'pushbutton', 'Position', [245, 10 100, 40], 'string', 'EdgeColor', 'callback', patchcolor)
    uicontrol(fi, 'Style', 'pushbutton', 'Position', [10, 10 100, 40], 'string', 'Apply', 'callback', applyfun)
end

function replotPatch(varargin)

    [d, fun, subset] = getData(); %#ok<ASGLU>
    set(0, 'CurrentFigure', mainFig)

    plotAsVector = any(ishandle(vecToggle)) && strcmpi(get(vecToggle, 'State'), 'on');
    plotAsLine = any(ishandle(linePlotToggle)) && strcmpi(get(linePlotToggle, 'State'), 'on');
    plotAsLine = plotAsLine || G.griddim == 1;
    if strcmpi(get(caxisToggle, 'State'), 'on')
        cx = caxis();
        yscale = ylim();
    end
    if strcmpi(get(gridToggle, 'State'), 'on')
        deleteHandle(gridOutline)
        gridOutline = plotGrid(G, 'FaceColor', 'none', 'EdgeAlpha', 0.05,...
            'EdgeColor', 1 - get(0, 'DefaultAxesColor'));
    else
        if any(ishandle(gridOutline))
            deleteHandle(gridOutline);
        end
    end

    if any(ishandle(vh))
        delete(vh);
    end
    if plotAsLine
        deleteHandle(gridOutline);
        
        if numel(ph) == 1 && strcmpi(get(ph, 'type'), 'line')
            color = get(ph, 'color');
            linewidth = get(ph, 'linewidth');
            style = get(ph, 'linestyle');
        else
            color = [0, 0, 0];
            linewidth = 2;
            style = '-';
        end
        deleteHandle(ph);
        useMarker = strcmpi(get(markerToggle, 'State'), 'on');
        if useMarker
            style = [style, 'o'];
        end
        if size(d, 2) > 1
            ph = plot(d, style, 'linewidth', linewidth);
        else
            ph = plot(d, style, 'linewidth', linewidth, 'color', color);
        end
    elseif plotAsVector
        if ishandle(ph)
            set(ph, 'Visible', 'off');
        end
        vh = plotCellVectorData(G, d, subset);
    else
        if numel(ph) <= 1 && ~isempty(ph) && ishandle(ph) && all(strcmpi(get(ph, 'Type'), 'Patch'))
            varg = {'EdgeAlpha', get(ph, 'EdgeAlpha'),...
                    'EdgeColor', get(ph, 'EdgeColor'),...
                    'FaceAlpha', get(ph, 'FaceAlpha'),...
                    'LineStyle', get(ph, 'LineStyle'),...
                    'LineWidth', get(ph, 'LineWidth')};
            if size(d, 1) == G.cells.num
                ph2 = plotCellData(G, d, subset, varg{:});
            else
                ph2 = plotNodeData(G, d, 'cells', find(subset), varg{:});
            end
            deleteHandle(ph);
            ph = ph2;
        else
            deleteHandle(ph);
            if size(d, 1) == G.cells.num
                if sum(subset)
                    ph = plotCellData(G, d, subset, 'EdgeColor', 'none', varargin{:});
                end
            else
                ph = plotNodeData(G, d, 'cells', find(subset), 'EdgeColor', 'none', varargin{:});
            end
        end
    end

    if strcmpi(get(zoomselToggle, 'State'), 'on')
        gc = G.cells.centroids(subset, :);
        mgc = min(gc);
        Mgc = max(gc);
        hgc = Mgc - mgc + sqrt(eps);
        axis([mgc(1) - .1*hgc(1)...
              Mgc(1) + .1*hgc(1)...
              mgc(2) - .1*hgc(2)...
              Mgc(2) + .1*hgc(2)...
              mgc(3) - .1*hgc(3)...
              Mgc(3) + .1*hgc(3)...
            ])
    end

    if strcmpi(get(caxisToggle, 'State'), 'on')
        caxis(cx);
        if plotAsLine
            % Interpret frozen colorbar as frozen yaxis
            ylim(yscale);
        end
    elseif min(size(d)) == 1
        if size(d, 1) == G.cells.num
            dsel = d(initialSelection);
        else
            dsel = d;
        end
        dsel = double(dsel(isfinite(dsel)));
        if any(dsel)
            lower = @(v) v - eps(v);
            upper = @(v) v + eps(v);
            caxis([lower(min(dsel)), upper(max(dsel))]);
        end
    end
    if opt.dynamicTitle
        if isempty(opt.time)
            s = sprintf('%s %d/%d', opt.title, ni, N);
        else
            s = sprintf('%s %d/%d: %s', opt.title,  ni, N, formatTimeRange(opt.time(ni), 2));
        end
    else
        s = opt.title;
    end
    if ~isempty(s)
        title(s);
    end
    cax = caxis();
    plotHistogram();
end

function [d, fun, subset] = getData()
    fun = @(x) x;
    if absDisplay
        fun = @(x) abs(fun(x));
    end
    if tenDisplay
        fun = @(x) 10.^fun(x);
    end
    if logDisplay
        fun = @(x) sign(fun(x)).*log10(abs(fun(x)));
    end
    subset = initialSelection;
    subset = subset & minmax.active;
    if ~isnan(slices)
        subset = subset & slices;
    end

    if ~isnan(histsubset)
        subset = subset & sum(histsubset, 2) > 0;
    end

    if ~isnan(logicalsubset.i)
        subset = subset & ismember(ijk{1}, logicalsubset.i);
    end

    if ~isnan(logicalsubset.j)
        subset = subset & ismember(ijk{2}, logicalsubset.j);
    end

    if ~isnan(logicalsubset.k)
        subset = subset & ismember(ijk{3}, logicalsubset.k);
    end

    if nonzeroDisplay
        tmp = sum(abs(data), 2);
        subset = subset & tmp > 1e-8*max(tmp);
    end

    if size(data, 2) == 3
        if ~any(ishandle(vecToggle))
            vecToggle = uitoggletool(ht, 'TooltipString', 'Plot as vector',...
                                         'ClickedCallback', @replotPatch, ...
                                         'cdata', geticon('vfield'));
        end
    elseif ishandle(vecToggle)
        delete(vecToggle);
    end
    d = fun(double(data));
    d = full(d);
end


end

function handleFigureResize(src, event)
    sel = findobj(src, 'Tag', 'mrst-datasetselector');

    oldunits_fig = get(src, 'Units');
    oldunits_sel = get(sel, 'Units');

    set(src, 'Units', 'pixels');
    set(sel, 'Units', 'pixels');

    figPos = get(src, 'Position');
    set(sel, 'Position', [0 0 figPos(3), 25]);

    set(src, 'Units', oldunits_fig);
    if iscell(oldunits_sel) && numel(oldunits_sel) > 1
        oldunits_sel = oldunits_sel{1};
    end
    set(sel, 'Units', oldunits_sel);
end

function cdata = geticon(name)
    this_dir = fileparts(mfilename('fullpath'));
    icon = fullfile(this_dir, 'icons', [name '.gif']);
    [cdata, map] = imread(icon);
    if islogical(cdata)
        cdata = uint8(255*cdata);
    end
    map((map(:,1)+map(:,2)+map(:,3)==3)) = NaN;
    cdata = ind2rgb(cdata,map);
end

function [s, e] = linkedSlider(fi, pos, md, Md, val, title)
    uicontrol(fi, 'Style', 'text', 'Position', [pos(1:2) 25, 25], 'string', title)
    size1 = [pos(1)+25, pos(2), 250, 25];
    sizee1 = [size1(1:2) + [size1(3) + 10, 0], [50 25]];
    cap = @(x) max(md, min(x, Md));

    e = uicontrol(fi, 'Style', 'edit', 'Value', val, 'Position', sizee1, 'String', sprintf('%.3f', val));
    fun = @(src, event) set(e, 'String', sprintf('%.3f', get(src, 'Value')));
    s = uicontrol(fi, 'Style', 'slider', 'Position', size1, 'Min', md, 'Max', Md, 'Value', val, 'Callback', fun);
    fun2 = @(src, event) set(s, 'Value', cap(sscanf(get(src, 'String'), '%f')));
    set(e, 'Callback', fun2);
end

function f = boolToOnOff(status)
    v = {'off', 'on'};
    f = v{status + 1};
end

function v = isPlottable(x)
    v = islogical(x) || isnumeric(x);
end
