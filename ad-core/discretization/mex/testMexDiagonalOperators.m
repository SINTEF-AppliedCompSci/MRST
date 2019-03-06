function results = testMexDiagonalOperators(model, varargin)
    N = model.operators.N;
    

    opt = struct('block_size', 5, 'iterations', 5, 'nc', max(N(:)));
    opt = merge_options(opt, varargin{:});
    
    results = opt;
    
    % Make test props
    nf = size(N, 1);
    nc = model.G.cells.num;
    d = rand(opt.nc, opt.block_size);
    flag = rand(nf, 1) > 0.5;
    cell_value = NewAD(rand(opt.nc, 1), DiagonalJacobian(d, size(d), []));
    % Perform tests
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %      Test face average      %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    f_avg = @(useMex) faceAverage(N, cell_value, useMex);
    [f_mex, f_matlab] = genFunctions(f_avg);
    [favg, results] = testFunction(f_mex, f_matlab, 'faceavg', 'Face average', opt, results);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Test single-point upwinding %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    upw = @(useMex) singlePointUpwind(flag, N, cell_value, useMex);
    [f_mex, f_matlab] = genFunctions(upw);
    [face_value, results] = testFunction(f_mex, f_matlab, 'upwind', 'Single-point upwind', opt, results);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %      Test gradient          %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    grad = @(useMex) twoPointGradient(N, cell_value, [], useMex);
    [f_mex, f_matlab] = genFunctions(grad);
    [gradient, results] = testFunction(f_mex, f_matlab, 'gradient', 'Two-point gradient', opt, results);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %      Prep divergence        %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tic();
    prelim = getMexDiscreteDivergenceJacPrecomputes(model);
    
    n1 = N(:, 1);
    n2 = N(:, 2);
    C  = sparse(N, [(1:nf)'; (1:nf)'], ones(nf,1)*[1 -1], nc, nf);
    gradMat = sparse((1:2*nf), [n1; n2], rldecode([-1; 1], [nf; nf]), 2*nf, nc);
    [~, sortedN] = sort(repmat(reshape(N, [], 1), 2, 1));
    I_base = [N(:, 1); N(:, 1); N(:, 2); N(:, 2)];

    sortIx = struct('C', C, 'J_sorted_index', sortedN, 'I_base', I_base);
    toc()
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %      Test divergence        %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    div = @(useMex) discreteDivergence([], N, face_value, nc, nf, sortIx, gradMat, prelim, useMex);
    [f_mex, f_matlab] = genFunctions(div);
    [divergence, results] = testFunction(f_mex, f_matlab, 'div', 'Discrete divergence', opt, results);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %      Test DivAcc            %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    f_mex = @() discreteDivergence(cell_value, N, face_value, nc, nf, sortIx, gradMat, prelim, true);
    f_matlab = @() cell_value + discreteDivergence([], N, face_value, nc, nf, sortIx, gradMat, prelim, false);
    [divacc, results] = testFunction(f_mex, f_matlab, 'accdiv', 'Accumulation + divergence', opt, results);
end

function [out, results] = testFunction(fn_mex, fn_mat, shortname, name, opt, results)
    its = opt.iterations;
    timer = tic();
    for i = 1:its
        matlab = fn_mat();
    end
    t_m = toc(timer);
    
    timer = tic();
    for i = 1:its
        mex = fn_mex();
    end
    t_c = toc(timer);
    out = matlab;
    v_error = norm(value(matlab) - value(mex))./norm(value(matlab));
    if ~issparse(mex.jac{1})
        mex.jac{1} = mex.jac{1}.sparse();
    end
    if ~issparse(matlab.jac{1})
        matlab.jac{1} = matlab.jac{1}.sparse();
    end
    j_error = norm(matlab.jac{1} - mex.jac{1}, inf)./norm(matlab.jac{1}, inf);
    fprintf('********************************************************\n');
    fprintf('* %s\n* Matlab: %1.2gs, MEX: %1.2gs (%1.2f speedup)\n', name, t_m/its, t_c/its, t_m/t_c);
    fprintf('Value error %1.2g, Jacobian error %1.2g\n', v_error, j_error)
    fprintf('********************************************************\n');
    
    if nargout > 1
        results.(shortname) = struct('t_matlab', t_m/its, 't_mex', t_c/its, 'description', name);
    end
end

function [f_mex, f_matlab] = genFunctions(fn)
    f_mex = @() fn(true);
    f_matlab = @() fn(false);
end