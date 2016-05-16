function A = polyArea3D(p)
% Computes the surface area of a polygon defined by the set of points 'p'
% in 3D

assert(size(p,1)>=3,'3 or more points needed to define a plane!');
assert(all(iscoplanar(p)),'All Points passed must be coplanar');

pdash = [p;p(1,:)];
diffp = diff(pdash,1);
normal = cross(diffp(1,:), diffp(2,:));
normal = normal/norm(normal);
A = 0;
for i = 1:size(p,1);
    A = A + cross(pdash(i,:),pdash(i+1,:));
end
A = 0.5*dot(normal,A);
return