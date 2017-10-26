function x = callAMGCL(A, b, varargin)
    opt = struct('coarsening',   'smoothed_aggregation', ...
                 'preconditioner', 'amg', ...
                 'relaxation',   'spai0', ....
                 'solver',       'bicgstab',...
                 'maxIterations', 0, ...
                 'tolerance',     1e-6);
    opt = merge_options(opt, varargin{:});
    switch(opt.preconditioner)
        case 'amg'
            precond = 1;
        case 'relaxation'
            precond = 2;
        case 'dummy'
            precond = 3;
        otherwise
            error('Unknown preconditioner option.')
    end
    
    switch(opt.coarsening)
        case 'smoothed_aggregation'
            coarse = 1;
        case 'ruge_stuben'
            coarse = 2;
        case 'aggregation'
            coarse = 3;
        case 'smoothed_aggr_emin'
            coarse = 4;
        otherwise
            error('Unknown coarsening option.')
    end

    switch(opt.relaxation)
        case 'spai0'
            relax = 1;
        case 'gauss_seidel'
            relax = 2;
        case 'ilu0'
            relax = 3;
        case 'iluk'
            relax = 4;
        case 'ilut'
            relax = 5;
        case 'damped_jacobi'
            relax = 6;
        case 'spai1'
            relax = 7;
        case 'chebyshev'
            relax = 8;
        otherwise
            error('Unknown relaxation option.')
    end
    
    switch(opt.solver)
        case 'bicgstab'
            solver = 1;
        case 'cg'
            solver = 2;
        case 'bicgstabl'
            solver = 3;
        case 'gmres'
            solver = 4;
        case 'lgmres'
            solver = 5;
        case 'fgmres'
            solver = 6;
        case 'idrs'
            solver = 7;
        otherwise
            error('Unknown solver option.')
    end
    % Note the transpose...
    x = amgcl_matlab(A', b, 1e-6, opt.maxIterations, coarse, relax, solver, precond);
end