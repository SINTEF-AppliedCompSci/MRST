function newtbl = duplicatefield(tbl, fdcell)
    
    oldfd = fdcell{1};
    fd1 = fdcell{2}{1};
    fd2 = fdcell{2}{2};

    [a, fds] = convertTableToArray(tbl, {oldfd});
    a = [a(:, 1), a];
    fds = fds(2 : end);
    fds = {fd1, fd2, fds{:}};
    
    newtbl = convertArrayToTable(a, fds);
    
end
