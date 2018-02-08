function [x, err] = callAMGCL_cpr(A, b, block_size, varargin)
    % Call AMGCL
    assert(block_size > 0);
    opt = struct('isTransposed',   false, ...
                 'cellMajorOrder', false, ...
                 'tolerance',     1e-6, ...
                 'maxIterations',  0);
    [opt, cl_opts] = merge_options(opt, varargin{:});
        
    amg_opt = getAMGCLMexStruct(cl_opts{:});
    amg_opt.block_size = block_size;
    
    assert(mod(size(A, 1), block_size) == 0);
    if opt.cellMajorOrder
        if ~opt.isTransposed
            A = A';
        end
    else
        n = size(A, 1);
        ordering = getCellMajorReordering(n/block_size, block_size, n);
        A = A(ordering, ordering)';
        b = b(ordering);
    end
    [x, err] = amgcl_matlab(A, b, amg_opt, opt.tolerance, opt.maxIterations, 2);
    
    if ~opt.cellMajorOrder
        x(ordering) = x;
    end
end