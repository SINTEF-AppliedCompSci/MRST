function [pts,removed] = surfaceSufCond3D(pts, CC, CR)
% Enforces the sufficient surface condition
%
% SYNOPSIS:
%   [pts, removed] = surfaceSufCond3D(pts, CC, CR)
%
% PARAMETERS:
%   pts     - A nx3 array of points. 
%   CC      - A mx3 array of Ball centers.
%   CR      - A mx1 array of Ball radii. 
%
% RETURNS:
%   pts     - A kx3 array of points (k<=n). No points in the retured pts is
%             inside any of the balls. E.g. If a point pts(i,:) from the 
%             input is inside any of the m balls it is removed. 
%   removed - A nx1 logical array. removed(i) is true if and only if input 
%             point pts(i,:) is innside at least one of the balls.
%
% EXAMPLE:
%   pts = rand(600,3);
%   CC  = [0,0,0;1,1,1];
%   CR  = [0.7;0.7];
%   pts = surfaceSufCond3D(pts, CC, CR);
%   figure()
%   plot3(pts(:,1), pts(:,2), pts(:,3),'.','markersize',15)
%
% SEE ALSO:
%   compositePebiGrid2D, pebi, clippedPebi3D, clipGrid.

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2016 Runar Lie Berge. See COPYRIGHT.TXT for details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}  


TOL = 10*eps;

nc = size(CC,1);
np = size(pts,1);
if np == 0
  removed = [];
  return
end

CRSqr = CR.^2;
removed = zeros(np,1);
for i = 1:nc
  distSqr = sum(bsxfun(@minus, CC(i,:), pts).^2,2);
  removed = removed + (distSqr<CRSqr(i)-TOL);
end
removed = logical(removed);
pts = pts(~removed,:);
end
