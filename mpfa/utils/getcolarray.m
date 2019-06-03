function vect = getcolarray(u, celltbl, cellcoltbl, cellfds, coltbl)
    ind = (1 : celltbl.num)';
    map = setupTableMapping(celltbl, cellcoltbl, cellfds);
    ind = map*ind;
    vect = sparse(ind, cellcoltbl.coldim, u, celltbl.num, coltbl.num);
    vect = full(vect);
end
