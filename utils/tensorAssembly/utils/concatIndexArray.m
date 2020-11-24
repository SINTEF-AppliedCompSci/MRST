function tbl = concatIndexArray(tbl1, tbl2, fdnames, varargin)
    
    opt = struct('checkUnique', true);
    opt = merge_options(opt, varargin{:}); 
    
    inds1 = tbl1.gets(fdnames);
    inds2 = tbl2.gets(fdnames);

    inds = [inds1; inds2];
    
    if opt.checkUnique
        indstest = unique(inds, 'rows');
        assert(size(indstest, 1) == size(inds, 1), 'there are repeated indices');
    end    
    
    tbl = IndexArray([], 'inds', inds, 'fdnames', fdnames);
    
end
