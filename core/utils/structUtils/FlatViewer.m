classdef FlatJsonViewer

    properties

        flatjson
        columnnames

    end

    methods

        function fjv  = FlatJsonViewer(flatjson, columnnames)

            if isa(flatjson, 'table')
                columnnames = flatjson.Properties.VariableNames;
                flatjson = table2cell(flatjson);
                skip_column_name_set = true;
            else
                skip_column_name_set = false;                
            end
            
            fjv.flatjson = flatjson;

            if (nargin < 2) & ~skip_column_name_set
                ncol = size(flatjson, 2);
                % default values for the column
                columnnames{1} = 'parameter name';
                if ncol == 2
                    columnnames{2} = 'parameter value';
                else
                    columnnames{2} = 'first set';
                    columnnames{3} = 'second set';
                    if ncol == 4
                        columnnames{4} = 'comparison';
                    end
                end
            end

            fjv.columnnames = columnnames;

        end

        function T = getTable(fjv)

            T = cell2table(fjv.flatjson, 'VariableNames', fjv.columnnames);
            
        end

        function T = print(fjv, varargin)

            opt = struct('filter', [], ...
                         'filename', [], ...
                         'print', true);

            opt = merge_options(opt, varargin{:});

            if ~isempty(opt.filter)
                fjv = fjv.filter(opt.filter);
            end

            T = cell2table(fjv.flatjson, 'VariableNames', fjv.columnnames);

            if opt.print
                display(T, 'view');
            end

            if ~isempty(opt.filename)
                writetable(T, opt.filename);
            end

        end


        function sortedfjv = sort(fjv, orderdesc)

            columnnames = fjv.columnnames;
            flatjson    = fjv.flatjson;

            ncol = numel(columnnames);

            if ischar(orderdesc)
                orderdesc = {orderdesc};
            end

            for iorder = numel(orderdesc) : -1 : 1

                r = regexprep(orderdesc{iorder}, ' +', '.*');
                ind = regexp(columnnames, r);
                ind = cellfun(@(res) ~isempty(res), ind);
                ind = find(ind);

                assert(numel(ind) == 1, 'regexp given is return too many match for column name');

                [~, ia] = sort(flatjson(:, ind));
                flatjson = flatjson(ia, :);

            end

            fjv.flatjson = flatjson;
            sortedfjv = fjv;

        end

        function filteredfjv = filter(fjv, filterdesc)

            if iscell(filterdesc{1})
                for ifilter = 1 : numel(filterdesc)
                    fjv = fjv.filter(filterdesc{ifilter});
                end
                filteredfjv = fjv;
                return
            end

            columnname  = filterdesc{1};
            filterval   = filterdesc{2};
            columnnames = fjv.columnnames;
            flatjson    = fjv.flatjson;

            r = regexprep(columnname, ' +', '.*');
            ind = regexp(columnnames, r);
            ind = cellfun(@(res) ~isempty(res), ind);
            ind = find(ind);
            assert(numel(ind) == 1, 'Ambiguous column name has been given');

            function res = strmatch(str)
                if regexp(str, r)
                    res = true;
                else
                    res = false;
                end
            end

            rowvals = flatjson(:, ind);

            if ischar(filterval)
                r = regexprep(filterval, ' +', '.*');
                filterval = @(str) strmatch(str);
            end

            ind = cellfun(@(res) filterval(res), rowvals);
            ind = find(ind);

            flatjson = flatjson(ind, :);

            fjv.flatjson = flatjson;
            filteredfjv = fjv;

        end

        function fjv = reorderColumns(fjv, frontColumNames)

            columnnames = fjv.columnnames;

            if iscell(frontColumNames)
                [isok, find] = ismember(frontColumNames, columnnames);
                assert(all(isok), 'some column names are not recognized.');
            else
                find = frontColumNames;
            end

            nind = true(numel(columnnames), 1);
            nind(find) = false;

            fjv.columnnames = horzcat(columnnames(find), ...
                                      columnnames(nind));

            fjv.flatjson = horzcat(fjv.flatjson(:, find), ...
                                   fjv.flatjson(:, nind));

        end
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
