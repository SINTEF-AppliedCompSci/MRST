function tbl = convertArrayToTable(A, fds)
    sz = size(A);
    nfds = numel(fds);
    
    assert(sz(2) == nfds, 'sizes are not matching');
    
    for ifield = 1 : nfds
        tbl.(fds{ifield}) = A(:, ifield);
    end
    tbl.num = sz(1);
end
