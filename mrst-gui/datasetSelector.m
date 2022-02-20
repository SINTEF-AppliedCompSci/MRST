function varargout = datasetSelector(G, datasets, varargin)
% Create dataset selector either as standalone or integrated in gui.
%
% SYNOPSIS:
%   % Select dataset and block until return
%   data = datasetSelector(G, datasets)
%
%   % Setup persistent selector with callback for integration with apps
%   datasetSelector(G, datasets, 'Parent', somehandle, 'Callback',
%   @callback)
%
% DESCRIPTION:
%  This function spawns a dataset selector which allows valid plotting
%  datasets to be selected. For instance, consider a state-object as
%  defined by initState:
%
%  The fields of this struct may be:
%  .s (Two columns of G.cells.num, suitable for plotting)
%  .pressure (One column of G.cells.num length, suitable for plotting)
%  .wellSol.pressure (Number of wells long column, not suitable for
%                     plotting)
%
%  datasetSelector deduces that only .s and pressure are suitable for
%  direct visualization, and presents these to the user. If multiple states
%  are provided as a struct array, these can be selected as well.
%
% REQUIRED PARAMETERS:
%   G    - Grid intended to visualization.
%
%   datasets - A valid dataset. Either a input suitable for plotCellData
%              directly, or a struct array of length 1 or more each
%              containing one or more fields which are valid input for
%              plotCellData, or have columns which are valid input for
%              plotCellData.
%
% OPTIONAL PARAMETERS:
%   Location  - Location supplied to parent gui for placement.
%
%   Parent    - Handle to gui component where the selector should be
%               embedded. If not provided (default), will spawn own window.
%
%   Setname   - Name to use for dataset. Defaults to the inputname of the
%               dataset, but this is not always meaningful with nested
%               calls.
%
%   Callback  - Function handle of the form @(data, N, name) where data is
%               input suitable for plotCellData, N is the index in the
%               datasets array and name is the name of the dataset.
%
%   Nofields  - Disables selecting / checking fields. 
%               Implies that both data and name will be [].
%   
%   Tag       - Tag to use for the panel containing the selector
%
%   pauseTime - Time (in seconds) between updates to plot when used for
%               dynamic playback. Note that this has no effect if the model
%               is too complicated for Matlab to render at this rate.
%
%
% RETURNS:
%   (OPTIONAL) Selected dataset.
%
% NOTE:
%   Behavior greatly dependent on number of output arguments. If output is
%   requested, the function will block for selection and close afterwards,
%   otherwise it will simply execute (possibly empty) callbback.
%
% SEE ALSO:
%   `plotToolbar`

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


opt = struct('Location', [0 0 1 1],...
             'Parent',   [], ...
             'Setname',  [], ...
             'active',   [], ...
             'N',        1, ...
             'Nofields', false,...
             'pauseTime',0.15, ...
             'Callback', [],...
             'Tag',      'mrst-datasetselector');

opt = merge_options(opt, varargin{:});

if isempty(opt.Parent)
    parent = figure('Toolbar','none', 'MenuBar', 'none');
    opt.Location = [0 0 1 .25];
else
    parent = opt.Parent;
end
ph = uipanel(parent, 'Position', opt.Location, ...
                     'Tag', opt.Tag);

if isnumeric(datasets)
    nsets = size(datasets, 2);
else
    if isstruct(datasets)
        nsets = numel(datasets);
        accessfcn = @(x) datasets(x);
    elseif iscell(datasets)
        nsets = numel(datasets);
        accessfcn = @(x) datasets{x};
    elseif strcmpi(class(datasets), 'ResultHandler')
        nsets = datasets.numelData();
        accessfcn = @(x) datasets{x};
    else
        error(['datasetSelector does not work for class ', class(datasets)])
    end
end
N = opt.N;

if isempty(opt.Setname) && ~ischar(opt.Setname)
    setname = inputname(2);
else
    setname = opt.Setname;
end

hasOutput = nargout > 0;
width = 1;
spos = 0;
if nsets > 1
    icons = {'<<', '<', 'Play', '>', '>>'};
    names = {'tostart', 'back', 'play', 'forward', 'toend'};
    increments = {-inf, -1, NaN, 1, inf};
    
    unused_width = 1.0;
    if (hasOutput)
        unused_width = unused_width - 0.2;
    end
    if (~opt.Nofields)
        unused_width = unused_width - 0.4;
    end
    h = unused_width/(numel(icons)+2);
    
    % Selection of iterable
    ith = uicontrol(ph, 'Style', 'popupmenu',...
                   'Units', 'normalized',...
                   'Value', N, ...
                   'String', arrayfun(@(i) sprintf('%d', i), 1:nsets, 'UniformOutput', false), ...
                   'Position', [0 0, 2*h, 1],...
                   'Callback', @selectionCallback...
                   );
    spos = spos + 2*h;
    width = width - 2*h;
    
    for i = 1:numel(icons)
        if isnan(increments{i})
            style = 'togglebutton';
        else
            style = 'pushbutton';
        end
        uicontrol(ph, 'Style', style,...
                   'Units', 'normalized',...
                   'String', icons{i},...
                   'Tag', [opt.Tag, '-', names{i}], ...
                   'Position', [spos 0, h, 1],...
                   'Callback', @(src, event) incrementDataset(src, event, increments{i})...
                   );
        spos = spos + h;
        width = width - h;
    end

end

if hasOutput
    uicontrol(ph, 'Style', 'pushbutton',...
                   'Units', 'normalized',...
                   'String', 'Apply', ...
                   'Position', [spos 0 .2 1],...
                   'Callback', @applySelection...
                   );
    spos = spos + .2;
    width = width - .2;
end

if (~opt.Nofields)
    cellfields = getStructFields(G, accessfcn(N), setname);
    assert(numel(cellfields) > 0, 'Number of plottable fields must be > 0');
    if isempty(opt.active) || numel(cellfields) < opt.active
        opt.active = 1;
    end
    structh = uicontrol(ph,  'Style', 'popupmenu',...
                                'Tag', 'mrst-activedataset', ...
                                'Value', opt.active, ...
                                'Units', 'normalized',...
                                'Position', [spos 0 max(width, eps) 1],...
                                'Callback', @applySelection, ...
                                'String', cellfields...
                               );
end

if nargout > 0
    uiwait(parent)
end

function selectionCallback(src, event)
    N = get(src, 'Value');
    if (~opt.Nofields)
        name = cellfields{get(structh, 'Value')};
        cellfields = getStructFields(G, accessfcn(N), setname);
        set(structh, 'String', cellfields)
        nameind = find(strcmp(cellfields, name));
        if isempty(nameind)
            nameind = 1;
        end
        set(structh, 'Value', nameind) 
    end
    applySelection();
end

function incrementDataset(src, event, step)
    if ~isnan(step)
        set(ith, 'Value', max(min(N+step, nsets), 1));
        selectionCallback(ith, []);
    else
        if N == nsets;
            N = 1;
        end
        for NN = N:nsets
            if ~get(src, 'Value')
                break;
            end
            set(ith, 'Value', NN);
            selectionCallback(ith, []);
            timer = tic();
            drawnow
            render_time = toc(timer);
            % We know the time it took to render the new model, try to fill
            % out the rest of the time by a explicit call to pause. If
            % rendering is slow, do nothing.
            pause_time = opt.pauseTime - render_time;
            if pause_time > 0
                pause(pause_time);
            end
        end
        set(src, 'Value', 0)
    end
end

function applySelection(src, event)
    if (opt.Nofields)
        name = [];
        d = [];
    else
        name = cellfields {get(structh, 'Value')};
        d = readStructField(accessfcn(N), name);
    end
    if hasOutput
        varargout{1} = d;
        varargout{2} = N;
        varargout{3} = name;
        close(parent)
    end
    if ~isempty(opt.Callback)
        opt.Callback(d, N, name)
    end
end

end
