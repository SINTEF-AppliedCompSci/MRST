function G = unitSquare(n)

[x, y] = meshgrid(0:1/(n-1):1, 0:1/(n-1):1);

size(x)

x(:,2:n-1) = x(:,2:n-1) + random('Normal', 0, 1/(3*(n-1)), n, n-2);
y(2:n-1,:) = y(2:n-1,:) + random('Normal', 0, 1/(3*(n-1)), n-2, n);

P = [x(:), y(:)];
t = delaunayn(P);

G = triangleGrid(P,t);
G = pebi(G);

