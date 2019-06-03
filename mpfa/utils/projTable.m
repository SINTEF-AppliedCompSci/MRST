function projtbl = projTable(tbl, fds)
    a = convertTableToArray(tbl, fds);
    nfds = numel(fds);
    a = a(:, (1 : nfds));
    a = unique(a, 'rows');
    projtbl = convertArrayToTable(a, fds);
end