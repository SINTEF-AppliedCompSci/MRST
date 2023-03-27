function G = makeTestPebiGrid()

N=7; M=5; [x,y] = ndgrid(0:N,0:M);
x(2:N ,2: M) = x(2:N,2:M) + 0.3*randn(N-1,M-1);
y(2:N ,2: M) = y(2:N,2:M) + 0.3*randn(N-1,M-1);
aG = pebi(triangleGrid([x(:) y(:)]));
G = makeLayeredGrid(aG, 3);
