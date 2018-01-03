function c = expandMatrixToCell(matrix, subset)
% Expand a matrix into cell arrays. Typical usage: Converting state
% representation of composition (as matrix) into AD-values (as cell array
% of columns vectors).
    if iscell(matrix)
        c = matrix;
        return
    end
    if nargin == 1
        subset = ':';
    end
    
    n = size(matrix, 2);
    c = cell(1, n);
    for i = 1:n
        c{i} = matrix(subset, i);
    end
end