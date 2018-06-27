function ordering = getCellMajorReordering(ncell, block_size, varargin)
    ncell_total = ncell*block_size;

    opt = struct('ndof', ncell_total, ...
                 'equation_ordering', [], ...
                 'cell_ordering', []);
    opt = merge_options(opt, varargin{:});

    ordering = (1:opt.ndof)';
    
    subs = 1:ncell_total;
    subs = reshape(subs, [], block_size)';
    if ~isempty(opt.cell_ordering)
        subs = subs(:, opt.cell_ordering);
    end
    if ~isempty(opt.equation_ordering)
        subs = subs(opt.equation_ordering, :);
    end
    subs = subs(:);
    ordering(1:ncell_total) = subs;
end