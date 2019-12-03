function [W] = createWellGridPoints3D(wellLines, rho)
% Creates a set of points folowing a well path.
%
% SYNOPSIS:
%   F = createWellPoints3D(wellLines, rho);
%
% PARAMETERS:
%   WellLines  - A nx1 cell array of well paths. Each element, 
%               wellLines{i}, must be a valid well path in 3D. The well  
%               path should be on the form [x1,y1,z1;x2,y2,z2;...] where 
%               each row gives a vertex of the well path:
%               .-----------------.----------------------.
%               (x1,y1,z1)        (x2,y2,z2)             (x3,y3,z3)
%   rho       - A function handle. rho{i}, is a function that sets the 
%               distance between the well points. If rho(x) = Const = 0.1
%               the distance between returned well points will be
%               approximately 0.1. 
%
%
% RETURNS:
%   W           - A struct with elements
%     W.pts     - A mx3 array of well points. The well points interpolates
%                 the given well paths, with a distance given by rho.
%     W.wellPos - A mapping from wellLines to W.nodes. 
%     W.nodes   - W.pts is not ordered randomly. Also two wells 
%                 might share one well point. To find the well points that
%                 belong to well i, we can use the following code:
%                 W.pts(W.nodes(wellPos(i):wellPos(i)-1),:) 
% 
% EXAMPLE:
%   wellLines = {[0,0,0;0.5,0,0;1,1,1],[0,0,0;0,1,0;1,1,1]};
%   rho = @(x) 0.1*ones(size(x,1),1);
%   W = createWellGridPoints3D(wellLines, rho);
%
%   figure()
%   hold on
%   wp1 = W.pts(W.nodes(W.wellPos(1):W.wellPos(2)-1),:);
%   wp2 = W.pts(W.nodes(W.wellPos(2):W.wellPos(3)-1),:);
%   plot3(wp1(:,1), wp1(:,2), wp1(:,3),'ro','markersize',5)
%   plot3(wp2(:,1), wp2(:,2), wp2(:,3),'b.','markersize',15)
%   plot3(wellLines{1}(:,1),wellLines{1}(:,2),wellLines{1}(:,3),'r')
%   plot3(wellLines{2}(:,1),wellLines{2}(:,2),wellLines{2}(:,3),'b')
%
% SEE ALSO:
%   createFaultGridPoints3D, clippedPebi3D, voronoi2mrst, 
%   createWellGridPoints

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2016 Runar Lie Berge. See COPYRIGHT.TXT for details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}  



W = zeros(0,3);
W.wellPos = ones(numel(wellLines)+1,1);
W.pts = [];
for i = 1:numel(wellLines)
  l = wellLines{i};
  
  lineDist = rho{i}((l(1,:) + l(end,:))/2);
  wellPts = interLinePath(l, rho{i}, lineDist, [0,0]);
  W.wellPos(i+1) = W.wellPos(i) + size(wellPts,1);
  W.pts = [W.pts; wellPts];
end

[W.pts,~,ic] = unique(round(W.pts*1e14)/1e14,'rows');
W.nodes = ic;

end