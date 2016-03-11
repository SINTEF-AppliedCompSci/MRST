function [Gp, bcp] = makePeriodicCartesianGrid(G)

% Uses the MRST module 'upscaling'
require upscaling

% Only make the dimensions with cartDims>1 periodic
dims  = find(G.cartDims>1);

bf    = getBoundaryFaces(G);
ndims = numel(dims);
bcl   = cell(1, ndims); % "left" faces
bcr   = cell(1, ndims); % "right" faces
for i = 1:ndims
    d = dims(i);
    bcl{i}.face = bf{d,1};
    bcr{i}.face = bf{d,2};
end

% Set pressure drop to zero for all directions as default
dp = cell(1, ndims);
dp(:) = {0};

% Create periodic grid
[Gp, bcp] = makePeriodicGridMulti3d(G, bcl, bcr, dp, 'dims', dims);
Gp.parent = G;

end

function bfaces = getBoundaryFaces(G)
% Finds the boundary faces of a grid.
% 
% This function finds the boundary faces of the given grid and returns
% the face indecies as a cell structure.

tags = [1 2; 3 4; 5 6];
boundaryFace = any(G.faces.neighbors == 0, 2);
ind = boundaryFace(G.cells.faces(:,1));
faceAndTag = G.cells.faces(ind, :);

bfaces = cell(G.griddim,2);
for i = 1:G.griddim
   bfaces{i,1} = faceAndTag(faceAndTag(:,2) == tags(i,1), 1);
   bfaces{i,2} = faceAndTag(faceAndTag(:,2) == tags(i,2), 1);
end

end

