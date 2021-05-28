function [results, models, info] = benchmarkAutoDiffBackends(model, backends, varargin)
%Benchmark routine that compares the execution speed and correctness of AD backends

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

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
    opt = struct('print', nargout == 0,...
                 'block_size', 5, ...
                 'iterations', 5, ...
                 'decimals',  3, ...
                 'verbose', mrstVerbose(), ...
                 'timingFn', @mean, ...
                 'walltime', true, ...
                 'testInterp', true, ...
                 'outputvalue', true, ...
                 'models', {{}}, ...
                 'timeit', false);
    opt = merge_options(opt, varargin{:});
    if opt.timeit
        opt.iterations = 1;
    end
    if ~iscell(backends)
        backends = {backends};
    end
    if iscell(model) && isstruct(model{1})
        % We are just getting the results for printing!
        results = model;
        info = struct();
        models = {};
    else
        dispif(opt.verbose, 'Benchmarking block size %d with %d iterations\n', opt.block_size, opt.iterations);
        [models, model] = gen_models(model, backends, opt);
        nb = numel(backends);
        results = cell(1, nb);
        if exist('rng', 'file')
            rng(0);
        end
        variables = cell(1, opt.block_size);
        [variables{:}] = deal(rand(model.G.cells.num, 1));

        for i = 1:nb
            dispif(opt.verbose, 'Benchmarking backend %d of %d...', i, nb);
            results{i} = perform_benchmark(models{i}, backends{i}, variables, opt);
            dispif(opt.verbose, ' ok.\n');
        end
        info = struct('iterations', opt.iterations, 'block_size', opt.block_size, ...
                       'nc', models{1}.G.cells.num, 'nf', size(models{1}.operators.N, 1));

    end
    
    results = compare_results(backends, results, opt);
end

function result = perform_benchmark(model, backend, var0, opt)
    op = model.operators;
    nf = size(op.N, 1);
    result = struct();
    
    var = var0;
    [var{:}] = backend.initVariablesAD(var0{:});
    x = var{1};
    if numel(var) > 1
        y = var{2};
    else
        y = 2*x;
    end
    if numel(var) > 2
        z = var{3};
    else
        z = 10*x;
    end
    for i = 2:numel(var)
        % Ensure all diagonals are present for first entry
        x = x.*var{i};
    end
    nc = numelValue(x);
     % ---------------- Cell -> cell operators -------------------------- %
    % Cell multiply
    fn = @() x.*y;
    result = bench('cell_xy', result, fn, opt);
    
    % Cell multiply with non-ad
    yv = value(y);
    fn = @() x.*yv;
    result = bench('cell_xv', result, fn, opt);
     
    % Cell multiply and add
    fn = @() x.*y + 2*z;
    [result, c_val] = bench('cell_xy_2z', result, fn, opt);
    
    % Diag mult
    fn = @() perform_diagMult(c_val, yv);
    result = bench('diagmult', result, fn, opt);
    
    % Diag product mult
    yv2 = yv;
    c_val2 = c_val;
    % Trigger copy-on-write
    yv2(1) = 3;
    c_val(1) = 3;
    fn = @() perform_diagProductMult(c_val, c_val2, yv, yv2);
    result = bench('diagproductmult', result, fn, opt);
    
    % 1D Interpolation
    if opt.testInterp
        n = [10, 50, 1000];
        for i = 1:numel(n)
            cv = value(c_val);
            X = linspace(min(cv) - 1, max(cv) + 1, n(i))';
            F = X.^2;
            if isprop(backend, 'useMex')
                if backend.useMex
                    fn = @() interpTableMEX(X, F, c_val);
                else
                    fn = @() interpTable(X, F, c_val);
                end
            else
                fn = @() interpTable(X, F, c_val);
            end
            result = bench(sprintf('interp1_%d', n(i)), result, fn, opt);
        end
    end
    
    
    fn = @() perform_diagMult(c_val, yv);
    result = bench('diagmult', result, fn, opt);
    
    % Sparse
    fn = @() perform_sparse(c_val);
    result = bench('sparse', result, fn, opt);
    
    % Subset (a few elements)
    subset_small = [1; ceil(nc/2); nc];
    % subs = struct('type', '()', 'subs', {{subset_small}});
    fn = @() perform_subset(c_val, subset_small);
    [result, v_small] = bench('subset_small', result, fn, opt);
    
    % Insert (a few elements)
    tmp = c_val;
    fn = @() perform_insertion(tmp, subset_small, v_small);
    result = bench('subasgn_small', result, fn, opt);

    % Subset (1% of elements)
    subset_med = 1:100:(nc-100);
    fn = @() perform_subset(c_val, subset_med);
    [result, v_med] = bench('subset_med', result, fn, opt);
    
    % Insert (1% of elements)
    tmp = c_val;
    fn = @() perform_insertion(tmp, subset_med, v_med);
    result = bench('subasgn_med', result, fn, opt);

    % Subset (many elements)
    subset_large = (1:5:nc)';
    fn = @() perform_subset(c_val, subset_large);
    [result, v_large] = bench('subset_large', result, fn, opt);
    
    % Insert (1% of elements)
    tmp = c_val;
    fn = @() perform_insertion(tmp, subset_large, v_large);
    result = bench('subasgn_large', result, fn, opt);
    
    % ---------------- Cell -> Face operators --------------------------- %
    % Face average
    fn = @() op.faceAvg(c_val);
    [result, f_val] = bench('faceavg', result, fn, opt);
    
    % Gradient
    fn = @() op.Grad(c_val);
    result = bench('Grad', result, fn, opt);

    % Face upstr
    flag = true(nf, 1);
    flag(2:2:end) = false;
    fn = @() op.faceUpstr(flag, c_val);
    result = bench('upw', result, fn, opt);


    % Face multiply
    fn = @() f_val.*f_val;
    result = bench('face_xy', result, fn, opt);
    
    % Face multiply with non-ad
    fv = value(f_val);
    fn = @() f_val.*fv;
    result = bench('face_xv', result, fn, opt);

    % Face multiply and add
    fn = @() f_val.*f_val + 2*f_val;
    result = bench('face_xy_2z', result, fn, opt);

    % ---------------- Face -> cell operators --------------------------- %
    % Div 
    fn = @() op.Div(f_val);
    result = bench('Div', result, fn, opt);
    
    % AccDiv
    fn = @() op.AccDiv(c_val, f_val);
    result = bench('AccDiv', result, fn, opt);
end

function [models, model] = gen_models(model, backends, opt)
    nb = numel(backends);
    models = opt.models;
    if isempty(models)
        dispif(opt.verbose, 'Building model...');
        % Need to set up
        if isnumeric(model)
            % We recieved dimensional inputs
            G = cartGrid(model);
            model = computeGeometry(G);
        end
        if isstruct(model)
            % We recieved a grid?
            G = model;
            rock = makeRock(G, 1, 1);
            model = ReservoirModel(G, rock, struct());
        end
        dispif(opt.verbose, ' ok.\n');
        models = cell(1, nb);
        for i = 1:nb
            dispif(opt.verbose, 'Setting up operators %d of %d...', i, nb);
            m = model;
            m.AutoDiffBackend = backends{i};
            m = m.validateModel();
            models{i} = m;
            dispif(opt.verbose, ' ok.\n');
        end
    else
        dispif(opt.verbose, 'Models already set up, verifying...');
        % Already set up, somehow
        nm = numel(models);
        assert(nm == nb);
        for i = 1:nb
            assert(isa(models{i}.AutoDiffBackend, class(backends{i})));
        end
        dispif(opt.verbose, ' ok.\n');
    end
    model = models{1};
end

function results = compare_results(backends, results, opt)
    reference = results{1};
    nb = numel(backends);
    totTime = zeros(1, nb);
    names = fieldnames(reference);
    maxlen = max(cellfun(@numel, names));
    dec = opt.decimals;
    for i = 1:nb
        if i == 1
            e = ' (baseline)';
        else
            e = '';
        end
        result = results{i};
        descr = backends{i}.getBackendDescription();
        fprintf('Backend #%d%s:\n%s:\n', i, e, descr);
        s = sprintf('%*s | %*s ', maxlen, 'Name', dec, 'Time (s)');
        if i > 1
            s = [s, sprintf('| %*s ', 5, 'Speedup')];
        end
        delim = repmat('-', 1, numel(s));
        disp(delim);
        disp(s);
        disp(delim);
        for j = 1:(numel(names)+1)
            last = j > numel(names);
            if last
                disp(delim);
                name = 'Total time';
                t_ref = totTime(1);
                t = totTime(i);
            else
                name = names{j};
                [ref, t_ref] = unpack(reference, name, opt);
                [new, t]     = unpack(result, name, opt);
            end
            fprintf('%*s | %*.*f', maxlen, name, 1, 6, t);
            if i > 1
                fprintf(' | %3.*f', dec-1, t_ref/t);
            end
            if opt.outputvalue && ~last
                e_val = compare_value(ref, new);
                e_jac = compare_jac(ref, new);
                if max(e_val, e_jac) > 1e-12
                    fprintf(' -> !!! Error %g in value, %g in Jacobian !!!', e_val, e_jac);
                end
            end
            fprintf('\n');
            if ~last
                totTime(i) = totTime(i) + t;
            end
        end
%         disp(repmat('-', 1, numel(s)));
%         fprintf('  -> Total time | %3.*fs', dec, totTime(i));
%         if i > 1
%             fprintf(' |%3.*f speedup)', dec - 1, totTime(1)/totTime(i));
%         end
        fprintf('\n');
    end
end

function [v, t] = unpack(x, name, opt)
    if isfield(x.(name), 'output')
        v = x.(name).output;
    else
        v = nan;
    end
    if opt.walltime
        t = x.(name).t_wall;
    else
        t = x.(name).t_cpu;
    end
    t = opt.timingFn(t);
end

function e = compare_value(ref, new)
    ref = value(ref);
    new = value(new);
    e = norm(ref - new, 2)/norm(ref);
end

function [ref, new] = get_jac(ref, new)
    both_jac = isa(ref.jac{1}, 'DiagonalJacobian') && ...
               isa(new.jac{1}, 'DiagonalJacobian');
    if both_jac && numel(ref.jac) == 1 && numel(new.jac) == 1
        ref = get_diag(ref);
        new = get_diag(new);
    else
        ref = get_sparse(ref);
        new = get_sparse(new);
    end
end

function e = compare_jac(ref, new)
    [ref, new] = get_jac(ref, new);

    e = norm(nonzeros(ref - new))/norm(nonzeros(ref));
end

function s = get_sparse(v)
    v = combineEquations(v);
    s = v.jac{1};
end

function s = get_diag(v)
    D = v.jac{1};
    s = D.diagonal;
    if D.rowMajor
        s = s';
    end
end


function [result, v] = bench(name, result, fn, opt)
    [t_wall, t_cpu, v] = time_op(fn, opt);
    if opt.outputvalue
        result.(name) = struct('t_wall', t_wall, 't_cpu', t_cpu, 'output', v);
    else
        result.(name) = struct('t_wall', t_wall, 't_cpu', t_cpu);
    end
end

function [wall_time, cpu_time, result] = time_op(fn, opt)
    n = opt.iterations;
    wall_time = zeros(1, n);
    cpu_time = zeros(1, n);
    for it = 1:n
        c = cputime();
        timer = tic();
        result = fn();
        wall_time(it) = toc(timer);
        cpu_time(it) = cputime() - c;
    end
end

function y = perform_subset(x, subset)
    y = x(subset);
end

function x = perform_insertion(x, subset, y)
    x(subset) = y;
end

function x = perform_sparse(x)
    for i = 1:numel(x.jac)
        if ~isnumeric(x.jac{i})
            x.jac{i}.sparse();
        end
    end
end

function x = perform_diagMult(x, v)
    D = [];
    for i = 1:numel(x.jac)
        [x.jac{i}, D] = diagMult(v, x.jac{i}, D);
    end
end

function x = perform_diagProductMult(x, y, v, w)
    D1 = [];
    D2 = [];
    for i = 1:numel(x.jac)
        [x.jac{i}, D1, D2] = diagProductMult(v, w, x.jac{i}, y.jac{i}, D1, D2);
    end
end
