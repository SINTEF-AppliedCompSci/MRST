function tbl = sortTable(tbl, fds);
    a = convertTableToArray(tbl, fds);
    a = sortrows(a);
    tbl = convertArrayToTable(a, fds);
end
