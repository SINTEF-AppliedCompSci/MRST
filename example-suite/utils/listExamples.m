function listExamples()
    % List examples in suite
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