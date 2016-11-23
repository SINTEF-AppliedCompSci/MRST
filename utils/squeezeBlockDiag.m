function A = squeezeBlockDiag(A, n, r, c)

[R, ~] = size(A);
N = numel(n);

if r == R
    A = A*repmat(speye(c),N,1);
else
    A = repmat(speye(r),1,N)*A;
end

end
