clc; clear all; close all;

load('basisGrid');

G = sortEdges(G);
G = computeVEM2DGeometry(G);

bEdg = find(any(G.faces.neighbors == 0,2));

n = size(P,1);
for i = 1:n

    tol = 1e-10;
    dofEdg1 = abs(G.faces.centroids(bEdg,1+l(i))          - ...
                 (G.faces.centroids(bEdg,2-l(i))*a(i) + b(i))) < tol;
    dofEdg2 = abs(G.faces.centroids(bEdg,1+l(mod(i,n)+1)) - ...
                 (G.faces.centroids(bEdg,2-l(mod(i,n)+1))*a(mod(i,n)+1) + b(mod(i,n)+1))) < tol;

    A = [P(mod(i,n)+1,:)-P(i,:);P(mod(i+1,n)+1,:)-P(i,:)]\[1;0];
    B = -P(i,:)*A;

    gD = @(X) X*A + B;

    bc = VEM2D_addBC([], G, bEdg(dofEdg1), 'pressure', gD);
    bc = VEM2D_addBC(bc, G, bEdg(dofEdg2), 'pressure', gD);
    bc = VEM2D_addBC(bc, G, bEdg(~dofEdg2 & ~dofEdg1), 'pressure', 0);

    sol = VEM2D(G, 0, 2, bc);

    figure;
    
    plotVEM2D(G, sol, 1, 'edgecolor', 'none');
    view(-37.5, 40)
    
end

figure;
fill(P(:,1), P(:,2), 'y', 'facealpha', .2)