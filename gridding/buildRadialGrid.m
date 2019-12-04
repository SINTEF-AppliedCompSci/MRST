function [G, t] = buildRadialGrid(p, nA, nR)
% Build the 2D radial grid from point and dimension definitions
%
% SYNOPSIS:
%   [G, t] = buildRadialGrid(p, nA, nR)
%
% PARAMETERS:
%  p     -   2D node coordinates, obey the logical numbering (angularly 
%            cycle fastest, then radially)
%  nA    -   Angular cell dimension
%  nR    -   Radial cell dimension
%
% RETURNS:
%  G - The 2D radial grid. The cells and nodes obey logical numbering 
%      (angularly cycle fastest, then radially). Each cell has four faces. 
%      The face types are:
%       Face 1:  Radial  -
%       Face 2:  Angular + 
%       Face 3:  Radial  +
%       Face 4:  Angular - 
%      If the points are generated from R+ to R-, the directions of face 1 
%      and 3 will be R+ and R-
%      If the points are generated in clockwise direction, the directions o
%      of face 2 and 4 will be A- and A+
%  t - Connectivity list
%
% EXAMPLE:
%   [nA, nR, rW, rM] = deal(40, 10, 2, 10);
%   th = linspace(0, 2*pi, nA+1); th = th(1:end-1);
%   r = logspace(log10(rW), log10(rM), nR+1);
%   [R, TH] = meshgrid(r, th);
%   [px, py] = pol2cart(TH(:), R(:));
%   p = [px(:), py(:)];
%   G = buildRadialGrid(p, nA, nR)
%   figure, axis equal, plotCellData(G, (1:G.cells.num)')
%   title('Cell array indices of the radial grid')
%  
% SEE ALSO:
%  `tessellationGrid` `pebi` `makeConnListFromMat`

    np = size(p,1);
    nd = reshape((1:np)', nA, nR+1);
    nd = [nd; nd(1,:)];
    t  = makeConnListFromMat(nd);
    G  = tessellationGrid(p, t);
    G.radDims = [nA, nR];
    [a, r] = ind2sub([nA, nR], (1:G.cells.num)');
    G.radIndices = [a, r];
    G.type = [G.type, { mfilename }];
end