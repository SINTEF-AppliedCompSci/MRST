classdef FlatStructViewer

    properties

        flatstruct
        columnnames

    end

    methods

        function fsv  = FlatStructViewer(flatstruct, columnnames)

            if isa(flatstruct, 'table')
                columnnames = flatstruct.Properties.VariableNames;
                flatstruct = table2cell(flatstruct);
                skip_column_name_set = true;
            else
                skip_column_name_set = false;                
            end
            
            fsv.flatstruct = flatstruct;

            if (nargin < 2) & ~skip_column_name_set
                ncol = size(flatstruct, 2);
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

            fsv.columnnames = columnnames;

        end

        function T = getTable(fsv)

            T = cell2table(fsv.flatstruct, 'VariableNames', fsv.columnnames);
            
        end

        function T = print(fsv, varargin)

            opt = struct('filter', [], ...
                         'filename', [], ...
                         'print', true);

            opt = merge_options(opt, varargin{:});

            if ~isempty(opt.filter)
                fsv = fsv.filter(opt.filter);
            end

            T = cell2table(fsv.flatstruct, 'VariableNames', fsv.columnnames);

            if opt.print
                display(T, 'view');
            end

            if ~isempty(opt.filename)
                writetable(T, opt.filename);
            end

        end


        function sortedfsv = sort(fsv, orderdesc)

            columnnames = fsv.columnnames;
            flatstruct    = fsv.flatstruct;

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

                [~, ia] = sort(flatstruct(:, ind));
                flatstruct = flatstruct(ia, :);

            end

            fsv.flatstruct = flatstruct;
            sortedfsv = fsv;

        end

        function filteredfsv = filter(fsv, filterdesc)

            if iscell(filterdesc{1})
                for ifilter = 1 : numel(filterdesc)
                    fsv = fsv.filter(filterdesc{ifilter});
                end
                filteredfsv = fsv;
                return
            end

            columnname  = filterdesc{1};
            filterval   = filterdesc{2};
            columnnames = fsv.columnnames;
            flatstruct    = fsv.flatstruct;

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

            rowvals = flatstruct(:, ind);

            if ischar(filterval)
                r = regexprep(filterval, ' +', '.*');
                filterval = @(str) strmatch(str);
            end

            ind = cellfun(@(res) filterval(res), rowvals);
            ind = find(ind);

            flatstruct = flatstruct(ind, :);

            fsv.flatstruct = flatstruct;
            filteredfsv = fsv;

        end

        function fsv = reorderColumns(fsv, frontColumNames)

            columnnames = fsv.columnnames;

            if iscell(frontColumNames)
                [isok, find] = ismember(frontColumNames, columnnames);
                assert(all(isok), 'some column names are not recognized.');
            else
                find = frontColumNames;
            end

            nind = true(numel(columnnames), 1);
            nind(find) = false;

            fsv.columnnames = horzcat(columnnames(find), ...
                                      columnnames(nind));

            fsv.flatstruct = horzcat(fsv.flatstruct(:, find), ...
                                   fsv.flatstruct(:, nind));

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
