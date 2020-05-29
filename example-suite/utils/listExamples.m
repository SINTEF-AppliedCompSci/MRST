function listExamples()
    % List examples in suite

%{
Copyright 2009-2020 SINTEF Digital, Mathematics & Cybernetics.

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

    pth = mrstPath('example-suite');
    [names, desc] = getExamplesInPath(fullfile(pth, 'examples-mrst'));
    
    printTable(names, desc, 'MRST');
    [names, desc] = getExamplesInPath(fullfile(pth, 'examples-user'));
    if isempty(names)
        names = {'user_example'};
        desc  = {'No user-defined examples'};
    end
    printTable(names, desc, 'USER-DEFINED');
end

%-------------------------------------------------------------------------%
function [names, descriptions] = getExamplesInPath(pth)
    list = dir(fullfile(pth, '*.m'));
    list = list(~[list.isdir]);
    [names, descriptions] = deal(cell(numel(list),1));
    for i = 1:numel(list)
        name = list(i).name(1:end-2);
        desc = feval(name);
        names{i} = name;
        descriptions{i} = desc;
    end
end

%-------------------------------------------------------------------------%
function printTable(names, descriptions, group)
    nn = max(cellfun(@numel, names));
    nd = max(cellfun(@numel, descriptions));
    
    fprintf('\n')
    nt =  2 + nn + 3 + nd + 2;
    bline = [repmat('=', 1, nt), '\n'];
    fprintf(bline);
    header = sprintf('%s EXAMPLES', group);
    nh = numel(header);
    fprintf('|');
    ntmp = round(nt/2 - nh/2);
    fprintf(repmat(' ', 1, ntmp));
    fprintf(header);
    fprintf(repmat(' ', 1, nt - ntmp - nh - 2));
    fprintf('|\n');
    fprintf(bline);
    for i = 1:numel(names)
        fprintf('| ')
        fprintf(['%', num2str(nn), 's'], names{i})
        fprintf(' | ');
        fprintf(['%-', num2str(nd), 's'], descriptions{i})
        fprintf(' |\n')
    end
    fprintf(bline);
    fprintf('\n');
end
