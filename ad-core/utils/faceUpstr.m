function xu = faceUpstr(flag, x, N, sz)
    if numel(flag) == 1
        flag = repmat(flag, size(N, 1), 1);
    end
    assert(numel(flag) == size(N, 1) && islogical(flag), ...
        'One logical upstream flag must'' be supplied per interface.');
    upCell       = N(:,2);
    upCell(flag) = N(flag,1);
    if isnumeric(x)
        % x is a simple matrix, we just extract values using the cells
        xu = x(upCell, :);
    else
        % x is likely AD, construct a matrix to achieve upstream weighting
        xu = sparse((1:sz(1))', upCell, 1, sz(1), sz(2))*x;
    end
end