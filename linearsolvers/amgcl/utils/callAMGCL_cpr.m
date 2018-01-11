function [x, err] = callAMGCL_cpr(A, b, block_size, varargin)
    % Call AMGCL
    assert(block_size > 0);
    opt = struct('coarsening',   'smoothed_aggregation', ...
                 'relaxation',   'spai0', ....
                 't_relaxation', 'spai0', ...
                 'solver',       'bicgstab',...
                 'isTransposed',  false, ...
                 'isCellOrdered', false, ...
                 'maxIterations', 0, ...
                 'tolerance',     1e-6);
    opt = merge_options(opt, varargin{:});
    
    coarse = translateOptionsAMGCL('coarsening', opt.coarsening);
    relax = translateOptionsAMGCL('relaxation', opt.relaxation);
    t_relax = translateOptionsAMGCL('relaxation', opt.t_relaxation);
    solver = translateOptionsAMGCL('solver', opt.solver);
    if ~opt.isTransposed
        A = A';
    end
    tic()
    [x, err] = amgcl_matlab_cpr(A, b, opt.tolerance, opt.maxIterations, coarse, relax, solver, t_relax, block_size);
    toc()
end