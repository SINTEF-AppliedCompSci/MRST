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
    names = cellfun(@(x) x.name,examples_mrst,'UniformOutput',false);
    desc = cellfun(@(x) x.description,examples_mrst,'UniformOutput',false);
    
    printTable(names, desc, 'MRST');
    
    examples_user = getExamplesInPath(fullfile(pth, 'examples-user'),'USER-DEFINED');
    
    if isempty(examples_user)
        names = {'user_example'};
        desc  = {'No user-defined examples'};
    else 
        names = cellfun(@(x) x.name,examples_user,'UniformOutput',false);
        desc = cellfun(@(x) x.description,examples_user,'UniformOutput',false);
    end
     
    printTable(names, desc, 'USER-DEFINED');
    examples = [examples_mrst; examples_user];
    
end

%-------------------------------------------------------------------------%
function examples = getExamplesInPath(pth,grp)
    list = dir(fullfile(pth, '*.m'));
    list = list(~[list.isdir]);
    examples = cell(numel(list),1);
    for i = 1:numel(list)
        name = list(i).name(1:end-2);
        desc = feval(name);
        examples{i}.name = name;
        examples{i}.description = desc;
        examples{i}.group = grp;
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
