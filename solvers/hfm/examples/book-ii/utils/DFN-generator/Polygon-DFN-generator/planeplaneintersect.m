function [point,direction]=planeplaneintersect(point1,normal1,point2,normal2)
% assumes the normal directions are not parallel

normal1=normal1/norm(normal1);
normal2=normal2/norm(normal2);

p1=dot(point1,normal1); % RHS of hessian form for plane 1
p2=dot(point2,normal2); % RHS of hessian form for plane 2

m=[normal1;normal2];
b=[p1;p2];

point=(m\b)';
direction=(null(m))'; direction=direction/norm(direction);

end