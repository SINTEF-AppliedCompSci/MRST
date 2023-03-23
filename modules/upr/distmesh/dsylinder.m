function d = dsylinder(p,x0, h,r)

d1    = @(p) sqrt((p(:,1)-x0(1)).^2 + (p(:,2)-x0(2)).^2) - r;
d2    = @(p) p(:,3) - x0(3);
d3    = @(p) -p(:,3) + (x0(3)-h);
d4    = @(p) sqrt(d1(p).^2 + d2(p).^2);
d5    = @(p) sqrt(d1(p).^2 + d3(p).^2);
da    = @(p) dintersect(dintersect(d1(p),d2(p)), d3(p));

d = zeros(size(p,1),1);

cond1 = d1(p)>0 & d2(p)>0;
cond2 = d1(p)>0 & d3(p)>0;
cond3 = ~cond1 & ~cond2;

d = zeros(size(p,1),1);
d(cond1) = d4(p(cond1,:));
d(cond2) = d5(p(cond2,:));
d(cond3) = da(p(cond3,:));
end