function F = removeSurfaceConflictSites3D(F)
% Removes conflicting surface sites from intersecting (or almost 
% intersecting) face constraints.
%
% SYNOPSIS:
%   Fr = removeSurfaceConflictSites3D(F)
%   
%
% PARAMETERS:
%   F         - Struct (as returned from surfaceSites3D) with 
%               elements:
%     F.f.pts   - Point coordinates.
%     F.f.Gs    - Grid spacing for each surface site. This is the distance
%                 between the two sites on oposite sides of the surface.
%     F.f.pri   - Priority of the surface sites.
%     F.f.l     - A mapping from surface sites to surfaces. Surface
%                 site F.f.pts(k,:) is generated from surface F.f.l(k).
%     F.c.CC    - The center coordinates of the balls used to create the
%                 sites F.f.pts.
%     F.c.CR    - The radii of the Balls used to create the sites F.f.pts.
%     F.l.f     - Map from surfaces to sites.
%     F.l.fPos  - Surface sites F.l.f(F.l.fPos(i):F.l.fPos(i+1)-1) belongs
%                 to surface i.
%     F.l.c     - Map from surfaces to F.c.CC.
%     F.l.cPos  - Ball F.l.c(F.l.cPos(i):F.l.cPos(i+1)-1) belongs to
%                 surface i
%
% RETURNS:
%   Fr        - Struct with same elements as F, but any conflict sites
%               have been removed. A site is a conflict site if it lies
%               inside any of the circles defined by (F.c.CC, F.c.R).
%               Circles belonging to surface i removes any conflict
%               points belonging to surface >= i+1.
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
%   F   = surfaceSites3D({fDt1,fDt2},{rho,rho});
%   Fh  = removeSurfaceConflictSites3D(F);
%   figure; hold on; axis equal
%   plot3(F.f.pts(:,1), F.f.pts(:,2), F.f.pts(:,3),'bo','markersize',5); 
%   plot3(Fh.f.pts(:,1), Fh.f.pts(:,2), Fh.f.pts(:,3),'r.','markersize',15);
%
% SEE ALSO:
%   compositePebiGrid3D, pebi, surfaceSites3D, 
%   lineSites3D.

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2016 Runar Lie Berge. See COPYRIGHT.TXT for details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}  

for i = 1:size(F.l.fPos)-2
  f1 = F.l.f(F.l.fPos(i):F.l.fPos(i+1)-1);
  c2 = F.l.c(F.l.cPos(i+1):F.l.cPos(i+2)-1);
  [~,removed] = surfaceSufCond3D(F.f.pts(f1,:),F.c.CC(c2,:),F.c.R(c2));
  F.f.pts(f1(removed),:) = []; 
  F.f.l(f1(removed)) = [];
  LIA = ismember(F.l.f,find(removed));
  sub = cumsum(LIA);
  F.l.f = F.l.f - sub;
  F.l.f(LIA) = [];
  F.l.fPos(2:end) = F.l.fPos(2:end) - sub(F.l.fPos(2:end)-1);
end
end
