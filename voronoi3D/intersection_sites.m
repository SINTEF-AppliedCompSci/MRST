function grid_1 = intersection_sites(intersections, ds)

grid_1 = {};
for i = 1:numel(intersections)
   int = intersections{i} ;
   p1 = int{1}(1,:);
   p2 = int{1}(2,:);
   v = p2 - p1;
   L = sqrt(sum(v.^2));
   v = v / L;
   n = round(L / ds);
   if L < 1.5*ds
       pts = zeros(0, 3);
   else
       t = linspace(ds/2, L - ds/2, n);
       steps = bsxfun(@times, v, t');
       pts = bsxfun(@plus, p1, steps);
   end  
   grid_1 = [grid_1, pts];
end


end