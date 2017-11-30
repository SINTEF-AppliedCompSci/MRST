function F = fixFaultIntersection(F)
% Removes conflicting fault points from intersecting (or almost 
% intersecting) faults.
%
% SYNOPSIS:
%   Fr = fixFaultIntersection(F)
%   
%
% PARAMETERS:
%   F         - Struct (as returned from createFaultGridPoints3D) with 
%               elements:
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
% RETURNS:
%   Fr        - Struct with same elements as F, but any conflict points
%               have been removed. A point is a conflict point if it lies
%               inside any of the circles defined by (F.c.CC, F.c.R).
%               Circles belonging to fault surface i removes any conflict
%               points belonging to fault surface >= i+1.
%
% EXAMPLE:  
%   nx = 10;
%   x = linspace(0,1,nx)';
%   z = x;
%   [X,Z] = meshgrid(x,z);
%   fDt1.ConnectivityList = delaunay([X(:),Z(:)]);
%   fDt2.ConnectivityList = delaunay([X(:),Z(:)]);
% 
%   fDt1.Points = [X(:), 0.5*X(:), Z(:)];
%   fDt2.Points = [X(:),-0.5*X(:), Z(:)];
%   rho = @(p) 1.5/nx*ones(size(p,1),1);
% 
%   F   = createFaultGridPoints3D({fDt1,fDt2},{rho,rho});
%   Fh  = fixFaultIntersection(F);
%   figure; hold on; axis equal
%   plot3(F.f.pts(:,1), F.f.pts(:,2), F.f.pts(:,3),'bo','markersize',5); 
%   plot3(Fh.f.pts(:,1), Fh.f.pts(:,2), Fh.f.pts(:,3),'r.','markersize',15);
%
% SEE ALSO:
%   compositePebiGrid3D, pebi, createFaultGridPoints3D, 
%   createWellGridPoints3D.

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2016 Runar Lie Berge. See COPYRIGHT.TXT for details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}  

for i = 1:size(F.l.fPos)-2
  f1 = F.l.f(F.l.fPos(i):F.l.fPos(i+1)-1);
  c2 = F.l.c(F.l.cPos(i+1):F.l.cPos(i+2)-1);
  [~,removed] = faultSufCond3D(F.f.pts(f1,:),F.c.CC(c2,:),F.c.R(c2));
  F.f.pts(f1(removed),:) = []; 
  F.f.l(f1(removed)) = [];
  LIA = ismember(F.l.f,find(removed));
  sub = cumsum(LIA);
  F.l.f = F.l.f - sub;
  F.l.f(LIA) = [];
  F.l.fPos(2:end) = F.l.fPos(2:end) - sub(F.l.fPos(2:end)-1);
end
end