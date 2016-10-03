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
% RETURNS:
%   W         - A mx3 array of well points. The well points interpolates
%               the given well paths, with a distance given by rho.
%
% EXAMPLE:
%   wellLines = {[0,0,0;0.5,0,0;1,1,1],[0,0,0;0,1,0;1,1,1]};
%   rho = @(x) 0.1*ones(size(x,1),1);
%   W = createWellGridPoints3D(wellLines, rho);
%
%   figure()
%   hold on
%   plot3(W(:,1), W(:,2), W(:,3),'.')
%   plot3(wellLines{1}(:,1),wellLines{1}(:,2),wellLines{1}(:,3),'r')
%   plot3(wellLines{2}(:,1),wellLines{2}(:,2),wellLines{2}(:,3),'b')
%
% SEE ALSO:
%   createFaultGridPoints3D, BallInt, clippedPebi3D, voronoi2mrst, 
%   createWellGridPoints

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2016 Runar Lie Berge. See COPYRIGHT.TXT for details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}  



W = zeros(0,3);
for i = 1:numel(wellLines)
  l = wellLines{i};
  lineDist = rho((l(1,:) + l(end,:))/2);
  W = [W; interLinePath(l, rho, lineDist, [0,0])];
end

end