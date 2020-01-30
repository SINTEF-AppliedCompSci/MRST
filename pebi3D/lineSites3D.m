function [W] = lineSites3D(cellConstraints, rho)
% Creates a set of points folowing a set of curves.
%
% SYNOPSIS:
%   F = lineSites3D(cellConstraints, rho);
%
% PARAMETERS:
%   WellLines  - A nx1 cell array of well paths. Each element, 
%               cellConstraints{i}, must be a valid path in 3D. The  
%               path should be on the form [x1,y1,z1;x2,y2,z2;...] where 
%               each row gives a vertex of the path:
%               .-----------------.----------------------.
%               (x1,y1,z1)        (x2,y2,z2)             (x3,y3,z3)
%   rho       - A function handle. rho{i}, is a function that sets the 
%               distance between the sites. If rho(x) = Const = 0.1
%               the distance between returned sites will be
%               approximately 0.1. 
%
%
% RETURNS:
%   W           - A struct with elements
%     W.pts     - A mx3 array of sites. The sites interpolates
%                 the given cell constraint, with a distance given by rho.
%     W.wellPos - A mapping from cellConstraints to W.nodes. 
%     W.nodes   - W.pts is ordered randomly. Also two paths 
%                 might share one well point. To find the sites that
%                 belong to cell constraint i, we can use the following code:
%                 W.pts(W.nodes(wellPos(i):wellPos(i)-1),:) 
% 
% EXAMPLE:
%   cellConstraints = {[0,0,0;0.5,0,0;1,1,1],[0,0,0;0,1,0;1,1,1]};
%   rho = @(x) 0.1*ones(size(x,1),1);
%   W = lineSites3D(cellConstraints, rho);
%
%   figure()
%   hold on
%   wp1 = W.pts(W.nodes(W.wellPos(1):W.wellPos(2)-1),:);
%   wp2 = W.pts(W.nodes(W.wellPos(2):W.wellPos(3)-1),:);
%   plot3(wp1(:,1), wp1(:,2), wp1(:,3),'ro','markersize',5)
%   plot3(wp2(:,1), wp2(:,2), wp2(:,3),'b.','markersize',15)
%   plot3(cellConstraints{1}(:,1),cellConstraints{1}(:,2),cellConstraints{1}(:,3),'r')
%   plot3(cellConstraints{2}(:,1),cellConstraints{2}(:,2),cellConstraints{2}(:,3),'b')
%
% SEE ALSO:
%   surfaceSites3D, clippedPebi3D, voronoi2mrstGrid3D, 
%   lineSites2D

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2016 Runar Lie Berge. See COPYRIGHT.TXT for details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}  



W = zeros(0,3);
W.wellPos = ones(numel(cellConstraints)+1,1);
W.pts = [];
for i = 1:numel(cellConstraints)
  l = cellConstraints{i};
  
  lineDist = rho{i}((l(1,:) + l(end,:))/2);
  wellPts = interLinePath(l, rho{i}, lineDist, [0,0], false);
  W.wellPos(i+1) = W.wellPos(i) + size(wellPts,1);
  W.pts = [W.pts; wellPts];
end

[W.pts,~,ic] = unique(round(W.pts*1e14)/1e14,'rows');
W.nodes = ic;

end
