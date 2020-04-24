function M = blockCSRtoSparse(A, verb)
% Convert a block CSR matrix to Matlab's scalar CSC matrix (block-wise
% ordering)
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