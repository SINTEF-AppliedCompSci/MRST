
a = [80 10 0;
     70 70 5;
     90 70 5;
     85 10 0];
 
b = rotatePlane(a, [0 0 1]);

p = mean(b);

diffp = diff(a,1);
normal = cross(diffp(1,:), diffp(2,:));

normal = normal/norm(normal);

diffp = diff(b,1);
normal2 = cross(diffp(1,:), diffp(2,:));

normal2 = normal2/norm(normal2);

pnew = rotatePlane([b;p],normal);