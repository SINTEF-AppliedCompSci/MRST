function varargout = editWells(G, W, rock, varargin)
%Edit well setup
%
% SYNOPSIS:
%   % Blocking call
%   W = editWells(G, W, rock)
%
% DESCRIPTION:
%   Interactive well editor.
%
% REQUIRED PARAMETERS:
%   G    - Grid
%
%   W    - Well structure as defined by for exampel addWell.
%
%   rock - Valid rock structure suitable for input to addWell.
%
% OPTIONAL PARAMETERS:
%   'callback'  - Callback function to execute once wells have changed.
%
%   edit  - Allow editing of wells (Note: Editing wells makes it easy to
%           produce ill-defined wells, so use with care.
%
%   plot  - Plot wells during the process.
%
%   add   - Allow adding of wells
%
%   figure - Figure where wells will be plotted.
%
%   parent - Parent figure for embedding of the panel.
%
%   remove - Allow removal of existing wells.
%
% RETURNS:
%   W - Updated well structure.
%
% NOTE:
%   This tool is meant for experienced users who know exactly what they
%   want to change. It at the moment quite easy to produce wells capable of
%   crashing most, if not all, the solvers in MRST.
%
% SEE ALSO:
%   `interactiveSelection`

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

hasOutput = nargout > 0;

opt = struct('edit', true, ...
            'parent', NaN, ...
            'plot', false, ...
            'figure', NaN, ...
            'add', true, ...
            'remove', true);

opt = merge_options(opt, varargin{:});

if ~ishandle(opt.parent)
    parent = figure('Toolbar','none', 'MenuBar', 'none');
else
    parent = opt.parent;
end

selectedWells = [];

if isempty(W)
    W = addWell([], G, rock, 1);
end

if isfield(W, 'lims')
    for wind = 1:numel(W)
        W(wind).lims = inf;
    end
end

wf = fields(W(1));
% Remove limits field as it is not handled, being a struct

isname = find(strcmpi(wf, 'name'));

ind = 1:numel(wf);
ind(isname) = 1;
ind(1) = isname;
wf = wf(ind);

columnname     = cellfun(@(x) lookupFieldName(x)   , wf, 'UniformOutput', false);
columnformat   = cellfun(@(x) lookupFieldType(W(1), x), wf, 'UniformOutput', false) .';
columneditable = repmat(opt.edit, 1, numel(columnname));

showToolbar = (opt.add || opt.remove);

if showToolbar
    bgh = uibuttongroup(parent, 'Units', 'normalized', 'Position', [0 .9 1 .1]);
    buttons = {'Add', 'Delete', 'Edit', 'Apply'};
    nb = numel(buttons);
    h = 1/nb;

    for ib = 1:nb
        uicontrol(bgh, ...
            'Units', 'normalized',...
            'String', buttons{ib},...
            'Callback', @buttonCallback, ...
            'Position', [h*(ib-1) 0 h 1])
    end
end

wellTable = uitable(parent,...
                  'Units','normalized',...
                  'Position', [0 0 1 1 - .1*showToolbar], ...
                  'Data', getData(W), ...
                  'ColumnName', columnname,...
                  'ColumnFormat', columnformat,...
                  'ColumnEditable', columneditable, ...
                  'CellSelectionCallback', @selectCallback,...
                  'CellEditCallback', @(src, event) updateWells() );

if hasOutput
    uiwait(parent)
    varargout{1} = W;
end

function selectCallback(src, event)
    selectedWells = event.Indices(:,1);
    plotCurrent()
end

function buttonCallback(src, event)
    switch get(src, 'String')
        case 'Delete'
            data = get(wellTable, 'Data');
            tmp = true(size(data, 1), 1);
            tmp(selectedWells, :) = false;
    
            %If all wells deleted, ignore, and warn user...
            if (sum(tmp) <= 1)
                warndlg('Please add at least two wells (one production, one injection)');
                return;
            end

            data = data(tmp, :);
            set(wellTable, 'Data', data);
            updateWells();
        case 'Add'
            cellInx = extractSubcellsInteractive(G, log10(rock.perm(:, 1)),...
                'selcount', 1, 'simple', true);
            if ~any(cellInx{1})
                return
            end
            W = addWell(W, G, rock, cellInx{1});
            set(wellTable, 'Data', getData(W));
        case 'Edit'
            if numel(selectedWells) > 1
                disp('Please select one well to edit...');
            end
            data = get(wellTable, 'Data');
            selection = extractSubcellsInteractive(G, log10(rock.perm(:, 1)), ...
                       'selectedCells', {[W(selectedWells).cells]}, ...
                       'selcount',      1);
            data{selectedWells, strcmpi(wf, 'cells')} = strjoin(selection{1});
            set(wellTable, 'Data', data);
            updateWells();
        case 'Apply'
            if hasOutput
                uiresume(parent);
                close(parent);
            end
            % Fire event...
    end
end

function descr = lookupFieldName(name)
    switch name
        case 'val'
            descr = 'Value';
        otherwise
            descr = name;
    end
end

function vtype = lookupFieldType(W, name)
    vtype = 'numeric';
    switch name
        case 'type'
            vtype = {'rate', 'bhp'};
        otherwise
            if isnumeric(W.(name)) && numel(W.(name) == 1)
                vtype = 'numeric';
            elseif ischar(W.(name))
                vtype = 'char';
            end
    end
end

function  data = getData(W)
    nw = numel(W);
    nf = numel(wf);
    data = cell(nw, nf);

    for i = 1:nw
        w = W(i);
        for j = 1:nf
            d = w.(wf{j});
            if isnumeric(d) || islogical(d)
                data{i, j} = strjoin(d);
            else
                data{i, j} = d;
            end
        end
    end
end

function updateWells()
    data = get(wellTable, 'Data');
    
    %If deleting well, clear the superfluous wells structures
    if (size(data, 1) < numel(W))
        indices = size(data, 1)+1:numel(W);
        W(indices) = [];
    end

    isval = strcmpi(wf, 'val');
    for i = 1:size(data, 1)
        d = cellfun(@(x) wellDataFromString(x), data(i, ind), 'UniformOutput', false);
        v = d{isval};
        if ischar(v)
            try
                v = eval(['[', d{isval}, ']']);
            catch
                v = NaN;
            end
            data{i, isval} = num2str(v);
        end
        W(i) = cell2struct(d .', wf(ind));
    end
    set(wellTable, 'Data', data);
    plotCurrent();
end

function plotCurrent()
    if ~opt.plot
        return
    end
    persistent f htop1 htext1 hs1 htop2 htext2 hs2
    if isempty(f)
        f = opt.figure;
    end
    old = [htop1 htext1 hs1 htop2 htext2 hs2];
    for i = 1:numel(old)
        if ishandle(old(i)); delete(old(i)); end;
    end

    if ~ishandle(f)
        f = figure;
    else
        set(0, 'CurrentFigure', f);
    end
    set(gca, 'ZDir', 'reverse')
    tmp = false(numel(W), 1);
    tmp(selectedWells) = true;

    if any(tmp)
        [htop1, htext1, hs1] = plotWell(G, W(tmp));
    end
    if ~all(tmp)
        [htop2, htext2, hs2] = plotWell(G, W(~tmp), 'Color', 'blue');
    end

end
end

function str = strjoin(iterable)
    si = size(iterable);
    assert(prod(si) == max(si));

    if si(1) > si(2)
        sep = '; ';
    else
        sep = ', ';
    end

    if ~iscell(iterable)
        iterable = arrayfun(@(x) num2str(x), iterable, 'UniformOutput', false);
    end
    if isempty(iterable)
        str = '';
    else
        str = [sprintf(['%s' sep],iterable{1:end-1}), iterable{end}];
    end
end

function data = wellDataFromString(data)
    if any(isletter(data))
        if ~all(strcmp(unique(data(isletter(data))), 'e'))
            return
        end
    end
    data = eval(['[', data, ']']);
end
