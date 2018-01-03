function [I, J, V, m, n] = getSparseArguments(M, ioffset, joffset)
% Get sparse matrix indices
    [I, J, V] = find(M);
    [m, n] = size(M);
    if nargin > 1
        I = I + ioffset;
        if nargin > 2
            J = J + joffset;
        end
    end
end