function [I, report, I_com] = cppMultiscaleBasis(CG, A, varargin)
% Create MsRSB basis using MEX code

% export LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libgomp.so.1

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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
    opt = struct('tolerance',  0.01, ...
                 'maxiter',    1000, ...
                 'omega',      2/3, ...
                 'basis',      [], ...
                 'verbose',    mrstVerbose,...
                 'writePath',  '', ...
                 'maxThreads', maxNumCompThreads(), ...
                 'doSolve',    true);
    opt = merge_options(opt, varargin{:});
    
    report = struct('t_setup', nan, 't_basis', nan);
    assert(CG.parent.cells.num == size(A, 1));
    if ~isfield(CG.cells, 'support_mex')
        error('Missing support_mex field in CG.cells. Did you call setupMexInteractionMapping?');
    end
    timer = tic();
    [offsets, support, types, I_com, fine, coarse] = getGridData(CG);
    report.t_setup = toc(timer);
    if ~isempty(opt.basis) && isfield(opt.basis, 'I_compressed')
        I_com = opt.basis.I_compressed;
    end
    
    if ~isempty(opt.writePath)
        % Error checking
        if ~exist(opt.writePath, 'dir')
            mkdir(opt.writePath)
        end
        
        inp = fullfile(opt.writePath, 'input');
        if ~exist(inp, 'dir')
            mkdir(inp);
        end
        [mat, j_index] = compressMatrix(A);
        % Write stuff
        writeGrid(CG, inp, offsets, support, types);
        writeMatrix(inp, j_index', mat');
        writeOperator(inp, I_com);
    end
    
    if opt.doSolve
        timer = tic();
        I_com = mex_cppMultiscaleBasis(offsets, support, types, A', I_com,...
                                   opt.tolerance/opt.omega, toIntegerValue(opt.maxiter), ...
                                   opt.omega, toIntegerValue(opt.maxThreads));
        report.t_basis = toc(timer);
    end
    I = sparse(fine, coarse, I_com, CG.parent.cells.num, CG.cells.num);
end

function [offsets, support, celltypes, I, fine, coarse] = getGridData(CG)
    sup = CG.cells.support_mex;
    offsets = sup.offsets;
    celltypes = sup.celltypes;
    support = sup.support;
    
    % Convert to 1-indexing
    ofs = double(offsets) + 1;
    I = zeros(size(support));
    for i = 1:CG.cells.num
        subs = ofs(i):(ofs(i+1)-1);
        I(subs) = CG.partition(support(subs, 1) + 1) == i;
    end
    
    fine = vertcat(sup.sorted_cell{:});
    coarse = rldecode((1:CG.cells.num)', cellfun(@numel, sup.sorted_cell));

end

function [mat, jj] = compressMatrix(A)
    n = size(A, 1);
    
    D = spdiags(1./diag(A), 0, n, n);
    [i_ix, j_ix, av] = find((D*A)');
    intx = i_ix ~= j_ix;
    i_ix = i_ix(intx);
    j_ix = j_ix(intx);
    av = av(intx);

    j2 = j_ix;

    pos = 1;
    t = tic();
    for i = 1:n
        v = j_ix(pos);
        ctr = 1;
        while v == i
            j2(pos) = ctr;
            ctr = ctr + 1;
            pos = pos + 1;
            if pos > numel(j_ix)
                break
            end
            v = j_ix(pos);
        end
    end
    
    mat = full(sparse(j_ix, j2, av))';
    jj = toIntegerIndex(full(sparse(j_ix, j2, i_ix)))';

    toc(t);
end

function v = toIntegerIndex(v)
    v = toIntegerValue(v - 1);
end

function v = toIntegerValue(v)
    v = int32(v);
end

function writeGrid(CG, fn, offsets, support, types)
    fh = fopen(fullfile(fn, 'info.txt'), 'w');
    fprintf(fh, '%d\r\n%d\r\n', CG.parent.cells.num, CG.cells.num);
    fprintf(fh, '%d ', offsets);
    fprintf(fh, '\r\n');
    fclose(fh);

    fh = fopen(fullfile(fn, 'support.txt'), 'w');
    fprintf(fh, '%d ', support');
    fprintf(fh, '\r\n');
    fclose(fh);

    fh = fopen(fullfile(fn, 'types.txt'), 'w');
    fprintf(fh, '%d ', types');
    fprintf(fh, '\r\n');
    fclose(fh);
end

function writeMatrix(fn, jj, mat)
    fh = fopen(fullfile(fn, 'sparsity.txt'), 'w');
    for i = 1:size(mat, 1)
        fprintf(fh, '%d ', jj(i, :));
        fprintf(fh, '\r\n');
    end
    fclose(fh);

    fh = fopen(fullfile(fn, 'matrix.txt'), 'w');
    fprintf(fh, '%d\r\n', size(mat, 2));
    for i = 1:size(mat, 1)
        fprintf(fh, '%1.8f ', mat(i, :));
        fprintf(fh, '\r\n');
    end
    fclose(fh);
end

function writeOperator(fn, I_init)
    fh = fopen(fullfile(fn, 'operator.txt'), 'w');
    fprintf(fh, '%1.8f ', I_init);
    fprintf(fh, '\r\n');
    fclose(fh); 
end