function [A, fds] = convertTableToArray(tbl, fdsfirst)
    
    fds = fieldnames(tbl);
    
    fdsToRemove = {fdsfirst{:}, 'ind', 'num'};
    for ifield = 1 : numel(fdsToRemove)
        fds = fds(~strcmp(fds, fdsToRemove{ifield}));
    end
    
    fds = {fdsfirst{:}, fds{:}};
    A = [];
    A = tbl.(fds{1});
    for ifield = 2 : numel(fds)
        A = [A, tbl.(fds{ifield})];
    end
     
end
