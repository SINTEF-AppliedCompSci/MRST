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
% OPTIONAL PARAMETERS (supplied in 'key'/value pairs ('pn'/pv ...)):
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
% RETURNS:
%   (OPTIONAL) Selected dataset.
%
% NOTE:
%   Behavior greatly dependent on number of output arguments. If output is
%   requested, the function will block for selection and close afterwards,
%   otherwise it will simply execute (possibly empty) callbback.
%
% SEE ALSO:
%   plotToolbar

%{
Copyright 2009-2014 SINTEF ICT, Applied Mathematics.

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
             'Callback', []);

opt = merge_options(opt, varargin{:});

if isempty(opt.Parent)
    parent = figure('Toolbar','none', 'MenuBar', 'none');
    opt.Location = [0 0 1 .25];
else
    parent = opt.Parent;
end
ph = uipanel(parent, 'Position', opt.Location, ...
                     'Tag', 'mrst-datasetselector');

if isnumeric(datasets)
    nsets = size(datasets, 2);
else
    nsets = numel(datasets);
end
N = 1;

if isempty(opt.Setname)
    setname = inputname(2);
else
    setname = opt.Setname;
end

width = 1;
spos = 0;
if nsets > 1
    % Selection of iterable
    ith = uicontrol(ph, 'Style', 'popupmenu',...
                   'Units', 'normalized',...
                   'Value', N, ...
                   'String', 1:nsets,...
                   'Position', [0 0, .1, 1],...
                   'Callback', @selectionCallback...
                   );
    icons = {'<<', '<', 'Play', '>', '>>'};
    increments = {-inf, -1, NaN, 1, inf};
    h = .3/numel(icons);
    for i = 1:numel(icons)
        if isnan(increments{i})
            style = 'togglebutton';
        else
            style = 'pushbutton';
        end
    uicontrol(ph, 'Style', style,...
                   'Units', 'normalized',...
                   'String', icons{i},...
                   'Position', [.10 + (i-1)*h 0, h, 1],...
                   'Callback', @(src, event) incrementDataset(src, event, increments{i})...
                   );
    end

    width = width - .4;
    spos = spos + .4;
end

hasOutput = nargout > 0;
if hasOutput
    uicontrol(ph, 'Style', 'pushbutton',...
                   'Units', 'normalized',...
                   'String', 'Apply', ...
                   'Position', [spos + width - .2 0 .2 1],...
                   'Callback', @applySelection...
                   );
    width = width - .2;
end

cellfields = getStructFields(G, datasets(N), setname);
structh = uicontrol(ph,  'Style', 'popupmenu',...
                            'Units', 'normalized',...
                            'Position', [spos 0 width 1],...
                            'Callback', @applySelection, ...
                            'String', cellfields...
                           );
if nargout > 0
    uiwait(parent)
end

function selectionCallback(src, event)
    name = cellfields{get(structh, 'Value')};
    N = get(src, 'Value');
    cellfields = getStructFields(G, datasets(N), setname);
    set(structh, 'String', cellfields)
    nameind = find(strcmpi(cellfields, name));
    if isempty(nameind)
        nameind = 1;
    end
    set(structh, 'Value', nameind)
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
        while N < nsets
            if ~get(src, 'Value')
                break;
            end
            set(ith, 'Value', N + 1);
            selectionCallback(ith, []);
            pause(.1)
        end
        set(src, 'Value', 0)
    end
end

function applySelection(src, event)
    name = cellfields {get(structh, 'Value')};
    d = readStructField(datasets(N), name);
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
