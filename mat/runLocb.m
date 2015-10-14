clc; clear all; close all;
clc; clear all; close all;

run('../../matlab/project-mechanics-fractures/mystartup.m')

nx = 20; ny = 20;
G = unitSquare(nx, ny);
G = mrstGridWithFullMappings(G);
G = computeGeometry(G);

g = @(X) -log(1./((X(:,1)+0.1).^2 + (X(:,2)+0.1).^2));
g = @(X) sin(X(:,1).*3*pi).*cos(X(:,2).*3*pi);

U = VEM2D(G,0,g);

fig1 = figure;
plotGridWithDofs(G);
fig2 = figure;
plotVEM(G, U, '');
fig3 = figure;
plotVEM(G, U, 'dof');
