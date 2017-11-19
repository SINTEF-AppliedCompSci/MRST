function plotDual_DFM(G, dual)
% Plot an implicitly defined dual grid for the multiscale finite volume method
%
% SYNOPSIS:
%   plotDual_DFM(G, dual)
%
% PARAMETERS:
%   G    - Grid structure
%   dual - Dual grid as defined by for example partitionUIdual.
%
%  Modified from plotDual.m to account for hybrid cells
%
%  Copyright 2013 IRIS AS
%
%  This file is licensed under the GNU General Public License v3.0.

RGB = @(r,g,b,i) [r/(i*255) g/(i*255) b/(i*255)];
plot(G.cells.centroids(dual.nn,1),G.cells.centroids(dual.nn,2),'*');
plotFractures(G, dual.ee);
plotGrid_DFM(G, 'facea',0,'edgea',.1);
axis tight
end
