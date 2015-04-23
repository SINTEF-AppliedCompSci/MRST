function D = formReconstructionMatrix(A, partition, keepDiag)
    if nargin == 2
        keepDiag = false;
    end
    N = size(A, 1);
    [i, j, d] = find(A);

    pi = partition(i);
    pj = partition(j);

    % Connections where both cells correspond to the same coarse block
    int = pi == pj;

    D  = sparse(i(int), j(int), d(int), N, N);
    if ~keepDiag
        D = D + diag(sum(A - D, 2));
    end
end
