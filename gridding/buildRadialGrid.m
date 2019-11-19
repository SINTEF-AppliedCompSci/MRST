function [G, t] = buildRadialGrid(p, nA, nR)
% Build the 2D radial grid from points and dimensions
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
%        faces(1): Radial- face
%        faces(3): Radial+ face
%        faces(2) and faces(4): Angular faces. The direction depends on the
%        cycle direction of input point p.
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
%  `tessellationGrid`

    np = size(p,1);
    nd = reshape((1:np)', nA, nR+1);
    nd = [nd; nd(1,:)];
    t  = makeConnListFromMat(nd);
    G  = tessellationGrid(p, t);
    G.radDims = [nA, nR];
    G.type = [G.type, { mfilename }];
end