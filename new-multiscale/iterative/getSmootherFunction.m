function fn = getSmootherFunction(varargin)
    opt = struct('type', 'jacobi', 'iterations', 5);
    opt = merge_options(opt, varargin{:});
    
    switch lower(opt.type)
        case 'jacobi'
            fn = @(A, b) getJacobiSmoother(A, b, opt.iterations);
        case 'ilu'
            fn = @(A, b) getILU0(A, b, opt.iterations);
        otherwise
            error('Unknown smoother')
    end
end

function smoother = getJacobiSmoother(A, ~, its)
    d = diag(A);
    Ar = A - diag(d);
    % D_inv = diag(1./d);
    n = size(A, 1);
    D_inv = spdiags(1./d, 0, n, n);
    % Apply a single pass of jacobi
    jac = @(d, x) D_inv*(d - Ar*x);
    smoother = @(d) loopfun(d, jac, its);
end

function smoother = getILU0(A, ~, its)
    setup = struct('type', 'nofill');
    
    [L, U] = ilu(A, setup);
    % Apply a single pass of jacobi
    jac = @(d, x) x + U\(L\(d - A*x));
    smoother = @(d) loopfun(d, jac, its);
end

function y = loopfun(d, smooth, its)
    y = zeros(size(d));
    for i = 1:its
        y = smooth(d, y);
    end
end