function flatjsonviewer = compareJsonStructs(jsonstruct1, jsonstruct2, name1, name2)
%%
% Compare two json structures, set a flag 'equal', 'different' or 'missing' on each entry. Returns the results as an instance of |FlatJsonViewer| for visualization.
%
    
    if nargin == 2
        name1 = 'jsonstruct1';
        name2 = 'jsonstruct2';
    end

    flatjsonviewer1 = flattenJsonStruct(jsonstruct1, 'doprint', false);
    flatjsonviewer2 = flattenJsonStruct(jsonstruct2, 'doprint', false);

    flatjson1 = flatjsonviewer1.flatjson;
    flatjson2 = flatjsonviewer2.flatjson;

    fd1 = flatjson1(:, 1);
    fd2 = flatjson2(:, 1);

    [fd, ia, ib] = intersect(fd1, fd2);

    flatjsondifferent = {};
    flatjsoncommon    = {};
    flatjsonmissing   = {};

    for ii = 1 : numel(ia)
        val1 = flatjson1{ia(ii), 2};
        val2 = flatjson2{ib(ii), 2};
        isequal = compareValue(val1, val2);
        entry = {fd{ii}, val1, val2};
        if isequal
            flatjsoncommon{end + 1} = entry;
        else
            flatjsondifferent{end + 1} = entry;
        end

    end

    ismissing1 = true(numel(fd2), 1);
    ismissing1(ib) = false;
    ismissing1 = find(ismissing1);
    for ii = 1 : numel(ismissing1)
        entry = {flatjson2{ismissing1(ii), 1}, NaN, flatjson2{ismissing1(ii), 2}};
        flatjsonmissing{end + 1} = entry;
    end

    ismissing2 = true(numel(fd1), 1);
    ismissing2(ia) = false;
    ismissing2 = find(ismissing2);
    for ii = 1 : numel(ismissing2)
        entry = {flatjson1{ismissing2(ii), 1}, flatjson1{ismissing2(ii), 2}, NaN};
        flatjsonmissing{end + 1} = entry;
    end

    flatjsoncommon    = vertcat(flatjsoncommon{:});
    flatjsondifferent = vertcat(flatjsondifferent{:});
    flatjsonmissing   = vertcat(flatjsonmissing{:});

    n                 = size(flatjsoncommon, 1);
    flatjsoncommon    = horzcat(flatjsoncommon, repmat({'equal'}, n, 1));
    n                 = size(flatjsondifferent, 1);
    flatjsondifferent = horzcat(flatjsondifferent, repmat({'different'}, n, 1));
    n                 = size(flatjsonmissing, 1);
    flatjsonmissing   = horzcat(flatjsonmissing, repmat({'missing'}, n, 1));

    flatjsonboth = vertcat(flatjsoncommon, flatjsondifferent, flatjsonmissing);

    flatjsonviewer = FlatJsonViewer(flatjsonboth, {'parameter name', name1, name2, 'comparison'});

end


function isequal = compareValue(val1, val2)

    if isstruct(val1) && isstruct(val2)
        % val1 and val2 are structs
        fds1 = fieldnames(val1);
        fds2 = fieldnames(val2);
        fds = intersect(fds1, fds2);
        if numel(fds) == numel(fds1) && numel(fds) == numel(fds2)
            % val1 and val2 have the same fields
            for ifd = 1 : numel(fds)
                % We compare the values of each field
                fd = fds{ids}
                subisequal = compareValue(val1.(fd), val2.(fd));
                if ~subisequal
                    isequal = false;
                    return
                end
            end
            isequal = true;
        else
            isequal = false;
        end
        return;
    end

    if iscell(val1) && iscell(val2)
        % val2 and val2 are cells
        if numel(val1) == numel(val2)
            % we compare the value of each cell
            for ival = 1 : numel(val1)
                subisequal = compareValue(val1{ival}, val2{ival});
                if ~subisequal
                    isequal = false;
                    return
                end
            end
            isequal = true;
        else
            isequal = false;
        end
        return
    end

    try
        isequal = eq(val1, val2);
    catch
        isequal = false;
    end

end


function printFunction_(cellarray, varargin)

    opt = struct('sortindex', [], ...
                 'columnorder', []);

    opt = merge_options(opt, varargin{:});

    nrow = size(cellarray, 1);
    ncol = size(cellarray, 2);

    strs = cell(nrow, ncol);
    ls   = zeros(ncol, 1);

    for ijson = 1 : nrow
        for col = 1 : ncol
            str = formattedDisplayText(cellarray{ijson, col});
            str = strtrim(str);
            str = regexprep(str, '\n', ',');
            str = convertStringsToChars(str);
            ls(col) = max(strlength(str), ls(col));
            strs{ijson, col} = str;
        end
    end

    if ~isempty(opt.sortindex)
        for isort = numel(opt.sortindex) : -1 : 1
            ind = opt.sortindex(isort);
            [~, ia] = sort(strs(:, ind));
            strs = strs(ia, :);
        end
    end

    if ~isempty(opt.columnorder)
        inds = (1 : ncol);
        ind1 = opt.columnorder;
        ind2 = setdiff(inds, ind1);
        inds = [ind1, ind2];
        strs = strs(:, inds);
    end

    formatstr = sprintf('%%-%ds', ls(1));
    for icol = 2 : ncol
        formatstr = sprintf('%s, %%-%ds', formatstr, ls(icol));
    end
    formatstr = sprintf('%s\\n', formatstr);

    for ijson = 1 : nrow
        str = strs(ijson, :);
        fprintf(formatstr, str{:});
    end

end



%{
Copyright 2021-2024 SINTEF Industry, Sustainable Energy Technology
and SINTEF Digital, Mathematics & Cybernetics.

This file is part of The Battery Modeling Toolbox BattMo

BattMo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

BattMo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with BattMo.  If not, see <http://www.gnu.org/licenses/>.
%}
