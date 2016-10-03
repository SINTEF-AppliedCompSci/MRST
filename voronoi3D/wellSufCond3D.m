function [pts,removed] = wellSufCond3D(pts, wellPts)
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
%   wellPts - A mx3 array of well points.
%
% RETURNS:
%   pts     - A kx3 array of points (k<=n). No points in the retured pts is
%             inside any of the balls. E.g. If a point pts(i,:) from the 
%             input is inside any of the m balls it is removed. 
%   removed - A nx1 logical array. removed(i) is true if and only if input 
%             point pts(i,:) is innside at least one of the balls.
%
% EXAMPLE:
% pts = rand(6000,3);
% x = linspace(0,1,5)';
% y = 0.5*ones(numel(x),1);
% z = 0.5*sin(x)+0.5;
% wellPts = [x,y,z];
% [~,removed] = wellSufCond3D(pts, wellPts);
% figure()
% hold on
% plot3(pts(removed,1), pts(removed,2), pts(removed,3),'.','markersize',15)
% plot3(wellPts(:,1),wellPts(:,2),wellPts(:,3),'.r');
% axis equal
%
% SEE ALSO:
%   compositePebiGrid, pebi, clippedPebi3D, clipGrid,
%   createWellGridPoints3D

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2016 Runar Lie Berge. See COPYRIGHT.TXT for details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}  
TOL = 10*eps;

removed = zeros(size(pts,1),1);
if size(wellPts,1)<=1  
  removed = logical(removed);  
  return
end
CC = (wellPts(1:end-1,:) + wellPts(2:end,:))/2;
RSqr = sum((wellPts(1:end-1,:) - wellPts(2:end,:)).^2,2)/4;

nc = size(CC,1);
np = size(pts,1);

for i = 1:nc
  distSqr = sum(bsxfun(@minus, CC(i,:), pts).^2,2);
  removed = removed + (distSqr<RSqr(i)-TOL);
end
removed = logical(removed);
pts = pts(~removed,:);
end






