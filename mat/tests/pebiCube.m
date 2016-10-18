function G = pebiCube(n,gridLim)
%   Generates a 3D pebi grid of n cells, covering the cube [0, gridLim(1)]
%   x [0, girdLim(2)] x [0, gridLim(3)].
%
%   SYNOPSIS:
%       G = voronoiCube(n, gridLim)
%
%   DESCRIPTION:
%       Generates 3D pebi grid with n cells covering the cube [0,
%       gridLim(1)] x [0, gridLim(2)] x [0, gridLim(3)] by definig the
%       convex boundary, and setting random seepoints inside the domain.
%       Based on [1]. Uses the function voronoi3D, see [1] and [2] for
%       details and copyright info.
%
%   REQUIRED PARAMETERS:
%       n       - Number of cells.
%       gridLim - Domain boundary.
%
%   RETURNS:
%       G   - MRST grid.
%
%   REFERENCES:
%       [1] - R. L. Berge. 'Unstructured pebi grids adapting to geological
%             feautres in subsurface reservoirs'. MA thesis. Norwegian
%             University of Science and Technology, 2016.
%       [2] - https://github.com/92runarlb/pebiGridding.git
%-----------------------------------------------------------------ØSK-2016-

% Sett the convex boundary
boundary = [0,          0         , 0         ; ...
            gridLim(1), 0         , 0         ; ...
            gridLim(1), gridLim(2), 0         ; ...
            0         , gridLim(2), 0         ; ...
            0         , 0         , gridLim(3); ...
            gridLim(1), 0         , gridLim(3); ...
            gridLim(1), gridLim(2), gridLim(3); ...
            0         , gridLim(2), gridLim(3)];
        
% Generate random seed points
pts = bsxfun(@times, rand(n,3), gridLim);

% Generate voronoi grid
G = voronoi3D(pts, boundary);

end

