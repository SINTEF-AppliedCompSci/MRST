function results = testMexDiagonalOperators(model, varargin)
%Undocumented Test Routine

%{
Copyright 2009-2019 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

    N = model.operators.N;
    G = model.G;

    opt = struct('print', nargout == 0, 'block_size', 5, 'iterations', 5, 'nc', max(N(:)), 'testSparse', true);
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
        d = cell_value.jac{1}.diagonal(:, i);
        ix = (1:nc)';
        J{i} = sparse(ix, ix, d, nc, nc);
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
    %       Test diagonal mult    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fn = @(x, y, z) x.*y + z;
    cell_mult = @(useMex) fn(cell_value, 2*cell_value, 3*cell_value);
    [f_mex, f_matlab] = genFunctions(cell_mult);
    f_sparse = @(varargin) fn(cell_value_sparse, 2*cell_value_sparse, 3*cell_value_sparse);
    [cellmultval, cellmultval_sparse, results] = testFunction(f_mex, f_matlab, f_sparse, 'cellmult', 'Multiply and add (cell)', opt, results);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  Test face diagonal mult    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cell_mult = @(useMex) fn(face_value, 2*face_value, 3*face_value);
    [f_mex, f_matlab] = genFunctions(cell_mult);
    f_sparse = @(varargin) fn(face_value_sparse, 2*face_value_sparse, 3*face_value_sparse);
    [facemultval, facemultval_sparse, results] = testFunction(f_mex, f_matlab, f_sparse, 'facemult', 'Multiply and add (face)', opt, results);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %      Prep divergence        %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if opt.print
        tic();
    end
    prelim = getMexDiscreteDivergenceJacPrecomputes(model);
    
    n1 = N(:, 1);
    n2 = N(:, 2);
    C  = sparse(N, [(1:nf)'; (1:nf)'], ones(nf,1)*[1 -1], nc, nf);
    gradMat = sparse((1:2*nf), [n1; n2], rldecode([-1; 1], [nf; nf]), 2*nf, nc);
    [~, sortedN] = sort(repmat(reshape(N, [], 1), 2, 1));
    I_base = [N(:, 1); N(:, 1); N(:, 2); N(:, 2)];

    sortIx = struct('C', C, 'J_sorted_index', sortedN, 'I_base', I_base);
    if opt.print
        toc()
    end
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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %    Test sparse()           %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    V_matlab = cell_value.jac{1};
    V_matlab.useMex = false;
    V_mex = V_matlab;
    V_mex.useMex = true;
    V_default = face_value_sparse.jac{1};
    f_mex = @() V_mex.sparse();
    f_matlab = @() V_matlab.sparse();
    f_sparse = @() V_default;
    [cellsparse, cellsparse_sparse, results] = testFunction(f_mex, f_matlab, f_sparse, 'sparse', 'Class -> Sparse', opt, results);
end

function [out, out_sparse, results] = testFunction(fn_mex, fn_mat, fn_sparse, shortname, name, opt, results)
    its = opt.iterations;
    [matlab, t_m] = perform_benchmark(fn_mat, its);
    [mex, t_c] = perform_benchmark(fn_mex, its);
    
    out = matlab;
    if issparse(matlab)
        v_error = 0;
    else
        v_error = norm(value(matlab) - value(mex))./norm(value(matlab));
    end

    if issparse(matlab)
        j_error = norm(matlab - mex, inf)./norm(matlab, inf);
    else
        if ~issparse(mex.jac{1})
            mex.jac{1} = mex.jac{1}.sparse();
        end
        if ~issparse(matlab.jac{1})
            matlab.jac{1} = matlab.jac{1}.sparse();
        end
        j_error = norm(matlab.jac{1} - mex.jac{1}, inf)./norm(matlab.jac{1}, inf);
    end
    
    
    if opt.testSparse
        [out_sparse, t_s] = perform_benchmark(fn_sparse, its);
    else
        out_sparse = nan;
        t_s = nan;
    end
    if opt.print
        fprintf('********************************************************\n');
        fprintf('* %s\n* Diagonal: %4fs, Diagonal-MEX: %4fs (%1.2f speedup)\n', name, t_m/its, t_c/its, t_m/t_c);
        if opt.testSparse
            fprintf('*   Sparse: %4fs, Diagonal-MEX: %4fs (%1.2f speedup)\n', t_s/its, t_c/its, t_s/t_c);
        end
        fprintf('Value error %1.2g, Jacobian error %1.2g\n', v_error, j_error)
        fprintf('********************************************************\n');
    end
    if nargout > 1
        results.(shortname) = struct('t_matlab', t_m/its, 't_mex', t_c/its, 't_sparse', t_s/its, 'description', name);
    end
end

function [out, t] = perform_benchmark(fn, its)
    timer = tic();
    bad = false;
    try
        for i = 1:its
            out = fn();
        end
    catch
        bad = true;
    end
    t = toc(timer);
    if bad
        t = nan;
        out = nan;
    end
end

function [f_mex, f_matlab] = genFunctions(fn)
    f_mex = @() fn(true);
    f_matlab = @() fn(false);
end
