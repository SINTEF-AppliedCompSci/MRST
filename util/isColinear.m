function [id] = isColinear(pts)
% isColinear tests if the points pts are colinear, that is, they lie on a
% straight line.
%
% SYNOPSIS:
%   id = isColinear(pts)
%
% PARAMETERS:
%   pts     A nx3 array of the set of points.
%
% RETURNS:
%   id      A logical that is true if the points pts lie on a straight line
%
% EXAMPLE:
%   pts1 = [0,0,0; 1,0,0;1,1,1];
%   pts2 = [0,0,0; 1,0,0;2,0,0];
%   isColinear(pts1)
%   isColinear(pts2)
%
% SEE ALSO:
%   calcNormals

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2016 Runar Lie Berge. See COPYRIGHT.TXT for details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}

id = rank(bsxfun(@minus,pts(2:end,:),pts(1,:)),1e-10)<2;
end