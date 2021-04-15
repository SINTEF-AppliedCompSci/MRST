function varargout = editStructGUI(data, varargin)
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

    if mod(numel(varargin), 2) == 1
        descr = varargin{1};
        varargin = varargin(2:end);
    else
        descr = [];
    end
    
    wantOutput = nargout > 0;
    yoffset = wantOutput*0.1;
    
    opt = struct('callback',    [], ...
        'figure',      [], ...
        'name',        '', ...
        'position',    [0, yoffset, 1, 1-yoffset]);
    opt = merge_options(opt, varargin{:});

    if isempty(opt.figure)
        opt.figure = figure('toolbar', 'none', 'name', opt.name);
    end


    tbl = uitable(opt.figure, 'units', 'normalized', ...
        'position',opt.position, ...
        'CellEditCallback',@editCallback);

    [d, lbl] = getData(data, descr);
    n_col = size(d, 2);
    edt = false(1, n_col);
    edt(2) = true;

    set(tbl, 'Data', d, 'ColumnName', lbl, 'ColumnEditable', edt);


    set(opt.figure, 'ResizeFcn', @(varargin) fixSizes(opt, n_col, tbl));


    function editCallback(hObject,callbackdata)
        try
            numval = eval(callbackdata.EditData);
        catch
            numval = nan;
        end
        r = callbackdata.Indices(1);
        c = callbackdata.Indices(2);
        tdata = get(hObject, 'Data');
        tdata{r,c} = numval;
        set(hObject, 'Data', tdata);

        data.(tdata{r, 1}) = numval;
        if ~isempty(opt.callback)
            % Send data as callback
            opt.callback(data);
        end
    end
    fixSizes(opt, n_col, tbl);
    
    if wantOutput
        uicontrol(opt.figure, 'units', 'normalized',...
                              'position', [0, 0, 1, yoffset], ...
                              'Style', 'pushbutton',...
                              'Callback', @(varargin) close(opt.figure), ...
                              'String', 'Done');
        
        uiwait(opt.figure);
        varargout{1} = data;
    end
end

function [d, lbl] = getData(data, descr)
    fnames = fieldnames(data);
    nvals = cellfun(@(x) numel(data.(x)), fnames);
    % Filter values with more than one element
    fnames = fnames(nvals == 1);

    fvals = cellfun(@(x) data.(x), fnames, 'UniformOutput', false);
    if isempty(descr)
        d = [fnames, fvals];
        lbl = {'Property', 'Value'};
    else
        n = numel(fnames);
        dsc = cell(n, 1);
        for i = 1:n
            fld = fnames{i};
            if isfield(descr, fld) && ischar(descr.(fld))
                dsc{i} = descr.(fld);
            else
                dsc{i} = '';
            end
        end

        d = [fnames, fvals, dsc];
        lbl = {'Property', 'Value', 'Description'};
    end
end

function fixSizes(opt, n_col, tbl)
    sz = get(opt.figure, 'Position');
    x = 0.91*sz(3);

    w = floor(x*1/n_col);

    width = cell(1, n_col);
    [width{:}] = deal(w);
    set(tbl, 'ColumnWidth', width);
end
