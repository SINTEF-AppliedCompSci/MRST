function results = testMexDiagonalOperators(model, varargin)
    N = model.operators.N;
    G = model.G;

    opt = struct('block_size', 5, 'iterations', 5, 'nc', max(N(:)), 'testSparse', true);
    opt = merge_options(opt, varargin{:});
    
    results = opt;
    
    % Make test props
    nf = size(N, 1);
    nc = model.G.cells.num;
    d = rand(opt.nc, opt.block_size);
    flag = rand(nf, 1) > 0.5;
    cell_value = GenericAD(rand(opt.nc, 1), DiagonalJacobian(d, size(d), []));
    
    % Make sparse version
    J = cell(1, opt.block_size);
    for i = 1:opt.block_size
        J{i} = spdiags(cell_value.jac{1}.diagonal(:, i), nc, nc);
    end
    ops_sparse = setupOperatorsTPFA(G, makeRock(G, 1, 1));
    cell_value_sparse = ADI(value(cell_value), J);
    % Perform tests
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %      Test face average      %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    f_avg = @(useMex) faceAverage(N, cell_value, useMex);
    [f_mex, f_matlab] = genFunctions(f_avg);
    f_sparse = @(varargin) ops_sparse.faceAvg(cell_value_sparse);
    
    [favg, favg_sparse, results] = testFunction(f_mex, f_matlab, f_sparse, 'faceavg', 'Face average', opt, results);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Test single-point upwinding %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    upw = @(useMex) singlePointUpwind(flag, N, cell_value, useMex);
    f_sparse = @(varargin) ops_sparse.faceUpstr(flag, cell_value_sparse);
    [f_mex, f_matlab] = genFunctions(upw);
    [face_value, face_value_sparse, results] = testFunction(f_mex, f_matlab, f_sparse, 'upwind', 'Single-point upwind', opt, results);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %      Test gradient          %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    grad = @(useMex) twoPointGradient(N, cell_value, [], useMex);
    [f_mex, f_matlab] = genFunctions(grad);
    f_sparse = @(varargin) ops_sparse.Grad(cell_value_sparse);
    [gradient, gradient_sparse, results] = testFunction(f_mex, f_matlab, f_sparse, 'gradient', 'Two-point gradient', opt, results);
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
    f_sparse = @(varargin) ops_sparse.Div(face_value_sparse);
    [divergence, div_sparse, results] = testFunction(f_mex, f_matlab, f_sparse, 'div', 'Discrete divergence', opt, results);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %      Test DivAcc            %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    f_mex = @() discreteDivergence(cell_value, N, face_value, nc, nf, sortIx, gradMat, prelim, true);
    f_matlab = @() cell_value + discreteDivergence([], N, face_value, nc, nf, sortIx, gradMat, prelim, false);
    f_sparse = @() ops_sparse.AccDiv(cell_value, face_value);
    [divacc, divacc_sparse, results] = testFunction(f_mex, f_matlab, f_sparse, 'accdiv', 'Accumulation + divergence', opt, results);
end

function [out, out_sparse, results] = testFunction(fn_mex, fn_mat, fn_sparse, shortname, name, opt, results)
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
    
    if opt.testSparse
        timer = tic();
        for i = 1:its
            out_sparse = fn_sparse();
        end
        t_s = toc();
    else
        out_sparse = nan;
        t_s = nan;
    end
    fprintf('********************************************************\n');
    fprintf('* %s\n* Diagonal: %1.2gs, Diagonal-MEX: %1.2gs (%1.2f speedup)\n', name, t_m/its, t_c/its, t_m/t_c);
    if opt.testSparse
        fprintf('*   Sparse: %1.2gs, Diagonal-MEX: %1.2gs (%1.2f speedup)\n', t_s/its, t_c/its, t_s/t_c);
    end
    fprintf('Value error %1.2g, Jacobian error %1.2g\n', v_error, j_error)
    fprintf('********************************************************\n');
    
    if nargout > 1
        results.(shortname) = struct('t_matlab', t_m/its, 't_mex', t_c/its, 't_sparse', t_s/its, 'description', name);
    end
end

function [f_mex, f_matlab] = genFunctions(fn)
    f_mex = @() fn(true);
    f_matlab = @() fn(false);
end