function A = sparseBlockDiag(A, n, d)
%   [A1 ]    [A1       ]
%   [...] -> [   ...   ]
%   [An ]    [       An]

nn = numel(n);
[r,c] = size(A);

assert((d == 1 && sum(n) == r) || ...
       (d == 2 && sum(n) == c)       );
   
if d == 1
    ii = repmat(1:r,1,c);
    jj = mcolon(1:c, (1:c)+nn*c-1, c); 
    jj = rldecode(jj, repmat(n,1,c), 2);
    A = sparse(ii,jj,A(:),r,nn*c);
    
else
    jj = repmat(1:c,1,r);
    ii = mcolon(1:r, (1:r)+nn*r-1, r);
    ii = rldecode(ii, repmat(n,1,r), 2);
    A = A';
    A = sparse(ii,jj,A(:),nn*r,c);
end

end