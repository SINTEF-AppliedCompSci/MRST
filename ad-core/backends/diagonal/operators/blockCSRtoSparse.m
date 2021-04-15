function M = blockCSRtoSparse(A, verb)
% Convert a block CSR matrix to Matlab's scalar CSC matrix (block-wise
% ordering)

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

    if nargin < 1
        verb = false;
    end
    [colNo, rowPtr, V, n] = deal(A.col_no, A.row_ptr, A.val, A.n);
    block_size = sqrt(size(V, 1));
    assert(floor(block_size) == block_size);
    assert(ceil(block_size) == block_size);

    rowNo = rldecode((0:(n-1))', diff(rowPtr)) + 1;
    colNo = double(colNo) + 1;
    
    ne = size(V, 2);
    A = repmat(1:block_size, block_size, 1);
    B = A';
    row = cell(1, ne);
    col = cell(1, ne);
    data = cell(1, ne);
    ofs = @(i) (i-1)*block_size;
    for i = 1:ne
        r = ofs(rowNo(i)) + A(:);
        c = ofs(colNo(i)) + B(:);
        v = V(:, i);
        row{i} = r;
        col{i} = c;
        data{i} = v;
        if verb
            fprintf('%d %d:\n', rowNo(i), colNo(i));
            disp(reshape(r, block_size, [])');
            disp(reshape(c, block_size, [])');
            disp(reshape(v, block_size, [])');
        end
    end

    Ir = vertcat(row{:});
    Jr = vertcat(col{:});
    Vr = vertcat(data{:});
    M = sparse(Ir, Jr, Vr, n*block_size, n*block_size);
end
