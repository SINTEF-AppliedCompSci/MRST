function [F] = createFaultGridPoints3D(faultTri, rho)
% Creates a set of points equidistant on both sides of a triangulation.
%
% The points are created by creating a ball centered at each vertex of the
% triangulation. For each triangle, we calculate the intersection of the 
% three associated balls. The equidistant points are placed at these
% intersections.
%
% SYNOPSIS:
%   F = createFaultGridPoints3D(faultTri, rho);
%
% PARAMETERS:
%   faultTri  - A nx1 cell array of triangulations. Each element, 
%               faultTri{i}, must be a valid triangulation of a surface in 
%               3D. The triangulation must contain the fields 
%               faultTri.Points and faultTri.ConnectivityList as described 
%               by the matlab function delaunayTriangulation.
%   rho       - A nx1 a cell array of functions. Each element, rho{i}, is a
%               function that sets the radius of the balls created around
%               each vertex in the triangulation.
%
% RETURNS:
%   F         - Struct with elements:
%     F.f.pts   - Point coordinates.
%     F.f.Gs    - Grid spacing for each fault point. This is the distance
%                 between the two points on oposite sides of the fault.
%     F.f.pri   - Priority of the fault points. The points belonging to
%                 faultTri{i} is given a priority i.
%     F.f.l     - A mapping from fault points to fault surfaces. Fault
%                 point F.f.pts(k,:) is generated from surface F.f.l(k).
%     F.c.CC    - The center coordinates of the balls used to create the
%                 points F.f.pts. This is equivalent to the verties of the
%                 triangulations in the cell arrray faultTri.
%     F.c.CR    - The radii of the Balls used to create the points F.f.pts.
%     F.l.f     - Map from triangulations, faultTri, to points.
%     F.l.fPos  - Fault points F.l.f(F.l.fPos(i):F.l.fPos(i+1)-1) belongs
%                 to triangulation faultTri{i}.
%     F.l.c     - Map from triangulations, faultTri, to F.c.CC.
%     F.l.cPos  - Ball F.l.c(F.l.cPos(i):F.l.cPos(i+1)-1) belongs to
%                 triangulation faultTri{i}.
%
% EXAMPLE:
%   [X,Y]     = meshgrid(0:0.2:1);
%   dt.Points = [X(:), Y(:)];
%   dt.ConnectivityList = delaunay([X(:), Y(:)]);
%   dt.Points = [dt.Points, zeros(size(dt.Points,1),1)];
%   rho = @(p) 0.3*ones(size(p,1),1);
%   F         = createFaultGridPoints3D({dt}, {rho});
%
%   figure()
%   hold on
%   patch('vertices',dt.Points,'faces',dt.ConnectivityList,'facealpha',0.3)
%   plot3(F.f.pts(:,1),F.f.pts(:,2),F.f.pts(:,3), '.','markersize',15)
%
% SEE ALSO:
%   createFaultGridPoints, BallInt, clippedPebi3D, voronoi2mrst

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2016 Runar Lie Berge. See COPYRIGHT.TXT for details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}  

% Initialize fault Structure
F.f.pts  = [];
F.f.Gs   = [];
F.f.pri  = [];
F.f.l    = [];
F.c.CC   = [];
F.c.R    = [];
F.l.fPos =  1;
F.l.cPos =  1;

for i = 1:numel(faultTri)
  dF        = faultTri{i};
  rFun      = rho{i};
  R         = rFun(dF.Points);

  CC        = reshape(dF.Points(dF.ConnectivityList',:)',9,[])';
  R         = reshape(R(dF.ConnectivityList',:),3,[])';
  [fPts,Gs] = ballInt(CC(:,1:3),R(:,1),CC(:,4:6),R(:,2),CC(:,7:9),R(:,3));
  p         = round(fPts*10^14)/10^14;
  [~, IA]   = unique(p,'rows');
  fPts      = fPts(IA,:);
  
  CC        = reshape(CC',3,[])';
  F.l.fPos  = [F.l.fPos; size(F.f.pts,1)+1+size(fPts,1)];
  F.l.cPos  = [F.l.cPos; size(F.c.CC,1)+1+size(CC,1)];
  F.f.pts   = [F.f.pts; fPts];
  F.c.CC    = [F.c.CC; CC];
  F.c.R     = [F.c.R; reshape(R',[],1)];
  F.f.pri   = [F.f.pri; i*ones(size(fPts,1),1)];
  F.f.Gs    = [F.f.Gs; Gs];
  F.f.l     = [F.f.l;  i*ones(size(fPts,1),1)];
end

F.l.f = (1:F.l.fPos(end)-1)';
F.l.c = (1:F.l.cPos(end)-1)';

end