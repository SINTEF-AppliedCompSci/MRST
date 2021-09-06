function examples = listExampleSuite()
    % List examples in suite

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

    pth = mrstPath('example-suite');
    
    examples_mrst = getExamplesInPath(fullfile(pth, 'examples-mrst'),'MRST');
    names = {examples_mrst.name};
    desc = {examples_mrst.description};
    
    printTable(names, desc, 'MRST');
    
    examples_user = getExamplesInPath(fullfile(pth, 'examples-user'),'USER-DEFINED');
    
    if isempty(examples_user)
        names = {'user_example'};
        desc  = {'No user-defined examples'};
    else 
        names = {examples_user.name};
        desc = {examples_user.description};
    end
     
    printTable(names, desc, 'USER-DEFINED');
    examples = [examples_mrst, examples_user];
    
end

%-------------------------------------------------------------------------%
function examples = getExamplesInPath(pth,grp)
    list = dir(fullfile(pth, '*.m'));
    list = list(~[list.isdir]);

    [~, name] = cellfun(@fileparts, {list.name}, 'UniformOutput', false);
    desc = cellfun(@feval, name, 'UniformOutput', false);
    examples = struct('name', name, 'description', desc, 'group', grp);
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
