function bfaces = getBoundaryFaces(G)
%Finds the boundary faces of a grid.
% 
% SYNOPSIS:
%   bfaces = getBoundaryFaces(G)
% 
% DESCRIPTION:
%   This function finds the boundary faces of the given grid and returns
%   the face indecies as a cell structure.
% 
% REQUIRED PARAMETERS:
%   G        - Grid structure.
% 
% RETURNS:
%   bfaces   - Boundary face indecies as a cell structure.
%                bfaces{d,1} - The face indecies for the first side of
%                              dimension d.
%                bfaces{d,2} - The face indecies for the second side of
%                              dimension d.
% 

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

