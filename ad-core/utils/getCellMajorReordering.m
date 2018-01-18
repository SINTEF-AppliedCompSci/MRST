function ordering = getCellMajorReordering(ncell, block_size, ndof)
    ncell_total = ncell*block_size;
    if nargin < 3
        ndof = ncell_total;
    end
    ordering = (1:ndof)';
    
    subs = 1:ncell_total;
    subs = reshape(subs, [], block_size)';
    subs = subs(:);
    ordering(1:ncell_total) = subs;
end