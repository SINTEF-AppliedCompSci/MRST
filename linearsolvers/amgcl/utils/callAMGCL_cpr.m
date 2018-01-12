function [x, err] = callAMGCL_cpr(A, b, block_size, varargin)
    % Call AMGCL
    assert(block_size > 0);
    opt = struct('coarsening',    'smoothed_aggregation', ...
                 'relaxation',    'spai0', ....
                 't_relaxation',  'spai0', ...
                 'solver',        'bicgstab',...
                 'isTransposed',   false, ...
                 'cellMajorOrder', false, ...
                 'maxIterations',  0, ...
                 'tolerance',      1e-6);
    opt = merge_options(opt, varargin{:});
    
    coarse = translateOptionsAMGCL('coarsening', opt.coarsening);
    relax = translateOptionsAMGCL('relaxation', opt.relaxation);
    t_relax = translateOptionsAMGCL('relaxation', opt.t_relaxation);
    solver = translateOptionsAMGCL('solver', opt.solver);
    
    
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
    [x, err] = amgcl_matlab_cpr(A, b, opt.tolerance, opt.maxIterations, coarse, relax, solver, t_relax, block_size);
    
    if ~opt.cellMajorOrder
        x(ordering) = x;
    end
end