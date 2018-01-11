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
        [ii, jj, vv] = find(A);

        subs = 1:n;
        subs = reshape(subs, [], block_size)';
        subs = subs(:);
        
        remap = zeros(n, 1);
        remap(subs) = 1:n;
        
        if opt.isTransposed
            A = sparse(remap(ii), remap(jj), vv, n, n);
        else
            A = sparse(remap(jj), remap(ii), vv, n, n);
        end
        
        b = b(subs);
    end
    tic()
    [x, err] = amgcl_matlab_cpr(A, b, opt.tolerance, opt.maxIterations, coarse, relax, solver, t_relax, block_size);
    toc()
    
    if ~opt.cellMajorOrder
        x = x(remap);
    end
end