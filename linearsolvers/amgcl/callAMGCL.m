function [x, err] = callAMGCL(A, b, varargin)
% Call AMGCL
    opt = struct('isTransposed',  false, ...
                 'tolerance',     1e-6, ...
                 'maxIterations',  0);

    [opt, cl_opts] = merge_options(opt, varargin{:});

    amg_opt = getAMGCLMexStruct(cl_opts{:});
    
    % Note the transpose...
    if ~opt.isTransposed
        A = A';
    end
    [x, err] = amgcl_matlab(A, b, amg_opt, opt.tolerance, opt.maxIterations, 1);
end