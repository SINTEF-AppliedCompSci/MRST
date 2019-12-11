function truth=circleshadowintersect(newcircle,oldcircle,shadowmult,tol)
% Checks whether new circle invades shadow zone of old circle. We assume
% that both circles are parallel. The shadowzone will be a cylinder around
% the old circle. The cylinder will have height of shadowmult*diameter and
% radius of (1+shadowmult)*radius.

h=oldcircle.radius*shadowmult;
r=oldcircle.radius*(1+shadowmult);

% First check that the normal distance between the circles exceeds 'h'. If
% this is exceeded, the new circle definitely does not enter the shadowzone
% of the old circle.
centnew=newcircle.center;
centold=oldcircle.center;
normalnew=newcircle.normal; normalnew=normalnew/norm(normalnew);
normalold=oldcircle.normal; normalold=normalold/norm(normalold);
normalave=0.5*(normalnew+normalold);
normaldist=abs(dot(centnew-centold,normalave));
if (normaldist-h)>tol
    truth=false; return;
end

% if normal distance is within h, then we check the center to center
% distance between the two circles. We project this distance onto a
% direction parallel to the circles and compare this to the sum of the new
% circle's radius and the shadowzone radius.
centdist=norm(centnew-centold);
projcentdist=sqrt((centdist^2)-(normaldist^2));

sumrad=newcircle.radius+r;

if (projcentdist-sumrad)>tol % if the projected distance is longer than the sum of radii, then the new circle does not intersect shadowzone
    truth=false; return;
else
    truth=true;return;
end





end