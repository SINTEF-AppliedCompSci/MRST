function M = blockCSRSparseToMatlabSparse(I, J, V, n, m)
    block_size = sqrt(size(V, 1));
    assert(floor(block_size) == block_size);
    assert(ceil(block_size) == block_size);
%     assert(size(V, 2) == J(end))
%     assert(size(V, 2) == numel(I))
%     assert(size(V, 1) == m)
%     assert(max(I) <= n - 1);
%     offset = 1;
%     ptr = double(I); 
%     j = double(J);

    rowNo = rldecode((0:(n-1))', diff(J)) + 1;
    colNo = double(I) + 1;
    
    if 1
        ne = size(V, 2);
        A = repmat(1:block_size, block_size, 1);
        B = A';
        
%         [B, A] = deal(A, B);
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
            if 0
                fprintf('%d %d:\n', rowNo(i), colNo(i));
                disp(reshape(r, block_size, [])');
                disp(reshape(c, block_size, [])');
                disp(reshape(v, block_size, [])');
            end
        end
    else
        row = cell(block_size, block_size);
        col = cell(block_size, block_size);

        for iB = 1:block_size
            for jB = 1:block_size
                row{iB, jB} = block_size*(rowNo-1) + iB;
                col{iB, jB} = block_size*(colNo-1) + jB;

    %             row{iB, jB} = rowNo + n*(iB - 1);
    %             col{iB, jB} = colNo + n*(jB - 1);
            end
        end
    end
    Ir = vertcat(row{:});
    Jr = vertcat(col{:});
    Vr = vertcat(data{:});
%     Vr = V(:);
    M = sparse(Ir, Jr, Vr, n*block_size, n*block_size);
end