function n = calcNormals(ID, pts)
% Returns the normals of the given triangles
%
% SYNOPSIS:
%   n = calcNormals(ID, pts)
%
% PARAMETERS:
%   ID      A mx3 array of triangle indices. Each row corresponds to one 
%           triangle. The vertices of triangle k is given by pts(ID(k,:)).
%           E.g., ID = convhull(pts)
%   pts     A nx3 array of the coordinats of the vertices
%
% RETURNS:
%   n       Normal of each triangle
%
% EXAMPLE:
%   pts = [0,0,0; 1,0,0;1,1,1];
%   n   = calcNormals([1,2,3],pts);
%
% SEE ALSO:
%   isColinear

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2016 Runar Lie Berge. See COPYRIGHT.TXT for details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}

tri = pts(ID',:);
x0  = tri(1:3:end,:);
x1  = tri(2:3:end,:);
x2  = tri(3:3:end,:);
n   = cross(x1 - x0, x2 - x0);
n   = bsxfun(@rdivide, n, sqrt(sum(n.^2,2)));
end