A = rand(2,10)
cols = [2,3,5];
As = sparseBlockDiag(A,cols,2);
full(As)
Ans = squeezeBlockDiag(As,cols, 2,36)
A-Ans

A = rand(10,2)
n = [2,3,5];
As = sparseBlockDiag(A,n,1);
full(As)
Ans = squeezeBlockDiag(As,n, 10,2)
A-Ans

A = rand(2,10)
n = [2,3,5];
As = sparseBlockDiag(A,n,2);
full(As)
Ans = squeezeBlockDiag(As,n, 2,10)
A-Ans

