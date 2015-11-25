clc;

run('../../matlab/project-mechanics-fractures/mystartup.m')


u = @(X) 2*sin(X(:,1)*pi).*cos(X(:,2)*pi) -log(1./((X(:,1)+0.1).^2 + (X(:,2)+0.1).^2));

G = cartGrid([10,10],[1,1]);
G = computeGeometry(G);
X = G.nodes.coords;
U = u(G.nodes.coords);

[x,y] = meshgrid(0:0.1:1, 0:0.1:1);
Z = griddata(X(:,1), X(:,2), U, x, y);
figure()
surf(x,y,Z)
xlabel('x')
max(Z)