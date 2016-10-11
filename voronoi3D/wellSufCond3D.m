function [pts,removed] = wellSufCond3D(pts, W)
% Enforces the sufficient well condition
%
% SYNOPSIS:
%   [pts, removed] = wellSufCond(pts)
%
% PARAMETERS:
%   pts     - A nx3 array of points that should be removed if they violate 
%             the sufficient well condition. A point pts(i,:) is violate
%             the sufficent well condition for well point wellPts(l,:) and
%             wellPts(k,:) if it is inside the ball centered at 
%             (wellPts(k,:) + wellpts(l,:))/2 with radius 
%             norm(wellPts(k,:) - wellpts(l,:))/2
%   W           - A struct with elements
%     W.pts     - A mx3 array of well points. The well points interpolates
%                 the given well paths, with a distance given by rho.
%     W.wellPos - A mapping from wellLines to W.nodes. 
%     W.nodes   - W.pts is not ordered randomly. Also two wells 
%                 might share one well point. To find the well points that
%                 belong to well i, we can use the following code:
%                 W.pts(W.nodes(wellPos(i):wellPos(i)-1),:) 
%   
% RETURNS:
%   pts     - A kx3 array of points (k<=n). No points in the retured pts is
%             inside any of the balls. E.g. If a point pts(i,:) from the 
%             input is inside any of the m balls it is removed. 
%   removed - A nx1 logical array. removed(i) is true if and only if input 
%             point pts(i,:) is innside at least one of the balls.
%
% EXAMPLE:
%   pts = rand(6000,3);
%   x = linspace(0,1,5)';
%   y = 0.5*ones(numel(x),1);
%   z = 0.5*sin(x)+0.5;
%   wellPts = {[x,y,z]};
%   W = createWellGridPoints3D(wellPts, @(x) 0.3*ones(size(x,1),1));
%   [~,removed] = wellSufCond3D(pts,W);
%   figure()
%   hold on
%   plot3(pts(removed,1), pts(removed,2), pts(removed,3),'.','markersize',15)
%   plot3(W.pts(:,1),W.pts(:,2),W.pts(:,3),'.r');
%   axis equal
%
% SEE ALSO:
%   clippedPebi3D, clipGrid, createWellGridPoints3D, faultSufCond3D

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2016 Runar Lie Berge. See COPYRIGHT.TXT for details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}  
TOL = 10*eps;

interWell = W.wellPos(2:end-1) - 1;
removed   = zeros(size(pts,1),1);
if size(W.pts,1)<=1  
  removed = logical(removed);  
  return
end
CC              = (W.pts(W.nodes(1:end-1),:) + W.pts(W.nodes(2:end),:))/2;
RSqr            = sum((W.pts(W.nodes(1:end-1),:) - W.pts(W.nodes(2:end),:)).^2,2)/4;
CC(interWell,:)   = [];
RSqr(interWell) = [];

nc = size(CC,1);

for i = 1:nc
  distSqr = sum(bsxfun(@minus, CC(i,:), pts).^2,2);
  removed = removed + (distSqr<RSqr(i)+TOL);
end
removed = logical(removed);
pts     = pts(~removed,:);
end