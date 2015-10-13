clc; clear all; close all;

run('../../matlab/project-mechanics-fractures/mystartup.m')

nx = 5; ny = 5;
G = unitSquare(nx, ny);
G = mrstGridWithFullMappings(G);
G = computeGeometry(G);

g = @(X) -log(1./((X(:,1)+0.1).^2 + (X(:,2)+0.1).^2));

U = VEM2D(G,0,g);

fig1 = figure;
plotGridWithDofs(G);
fig2 = figure;
plotVEM(G, U, '');

%   Make function for computing baricenters of all cells
%   Make function for detecting boundary dofs
%   Change VEM2D to accept force term f