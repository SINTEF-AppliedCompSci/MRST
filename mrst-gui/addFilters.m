function filter = addFilters(G, filters, datasets, callback)
%Filter cells based on dataset.
%
% SYNOPSIS:
%   addFilters(G, [], datasets, @(selection, filters) function)
%
% DESCRIPTION:
%
% REQUIRED PARAMETERS:
%
%   G - Grid data structure.
%
%   filters - Filters emenating from previous call to addFilters.
%
%   datasets - Dataset which is to be filtered.
%
%   callback - Function handle accepting (active, filters) interface.
%   Active is a logical mask suitable for plotting based on the filters and
%   filters a datastructure which allows filters to be reapplied in the
%   same gui.
%
%
% RETURNS:
%   Nothing
% NOTE:
%   Internal function for use by plotToolbar. Subject to change.
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

    parent = figure('Toolbar','none', 'MenuBar', 'none');

    if isempty(filters)
        filters = getFilters(getStructFields(G, datasets(1), inputname(3)));
    end
    columnname =   {'Dataset full name', 'Minimum value', 'Maximum value', 'Filter status'};
    columnformat = {'char', 'numeric', 'numeric', {'On', 'Off'}};
    columneditable =  [false true true true];

    minmaxh = uitable(parent,...
                      'Units','normalized',...
                      'Position', [0 0 1 1], ...
                      'Data', filters, ...
                      'ColumnName', columnname,...
                      'ColumnFormat', columnformat,...
                      'ColumnEditable', columneditable, ...
                      'ColumnWidth', 'auto', ...
                      'CellEditCallback', @applyFilters);

    function filters = getFilters(sfields)
        nf = numel(sfields);
        filters = cell(nf, 4);
        for i = 1:nf
            s = sfields{i};
            data = readStructField(datasets(1), s);
            if numel(data) > G.cells.num
                data = sqrt(sum(data.^2, 2));
            end
            filters{i, 1} = s;
            filters{i, 2} = min(data);
            filters{i, 3} = max(data);
            filters{i, 4} = 'On';
        end
    end

    function applyFilters(src, event)
        active = true(G.cells.num, 1);
        filters = get(minmaxh, 'data');

        for i = 1:size(filters, 1);
            if strcmpi(filters{i, 4}, 'off')
                continue;
            end
            data = readStructField(datasets(1), filters{i, 1});
            if size(data, 2) > 1
                data = sqrt(sum(data.^2, 2));
            end
            if isfield(G, 'nodes') && size(data, 1) == G.nodes.num
                if ~isfield(G.cells, 'nodes')
                    continue
                end
                cellNo = rldecode(1:G.cells.num, diff(G.cells.nodePos), 2) .';
                % Take average of node data to get rough approximation of
                % filter
                data = accumarray(cellNo, data(G.cells.nodes(:, 1)), [], @mean);
            end
            active = active & data(:,1) >= filters{i, 2} & data(:,1) <= filters{i, 3};
        end
        if ~isempty(callback)
            callback(active, filters);
        end
    end
    filter = @(varargin) applyFilters(varargin);
end
