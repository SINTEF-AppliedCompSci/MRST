function a = makexy(xcur)
for i =1:1:4
    newpos(i,:) = xcur(i,:) - xcur(1,:);
end
basis1 = (newpos(2,:))/norm(newpos(2,:)); % Gram Schmidt orthonormalisation assuming linear independence
x2d = newpos(3,:) - basis1*dot(newpos(3,:),basis1);
basis2 = x2d/norm(x2d);
basis3 = cross(basis1,basis2);
for i=1:1:4
    a(i,:) = [dot(newpos(i,:),basis1),dot(newpos(i,:),basis2),dot(newpos(i,:),basis3)];
end
d1 = -xcur(1,:);
d1d = [dot(d1,basis1),dot(d1,basis2),dot(d1,basis3)];
check = norm(d1) - norm(d1d) % norms must be equal 
