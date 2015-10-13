function G = unitDisc(n)

r = 0;
for rho = 0:0.1:1
    rTmp = rho.*ones(n*ceil(3.85*rho^2),1)
    r = [r;rTmp];
end
m = numel(r);

theta = randperm(m).*2.*pi./m;


P = [r.*cos(theta'), r.*sin(theta')];
t = delaunayn(P);

theta = 0:0.01:2*pi;
G = triangleGrid(P,t);
G = pebi(G);
hold on;
plotGrid(G)
plot(cos(theta), sin(theta))
axis equal
hold off


