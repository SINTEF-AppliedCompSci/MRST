function intersections = fault_intersections(faults)

intersections = {};
for f1 = 1:numel(faults)
    % first find intersection between fault i and all other faults
    for f2 = f1+1:numel(faults)
       int_pts = polygon_intersection(faults{f1}, faults{f2});
       if size(int_pts, 1) > 0
           intersections= {intersections{:}, {int_pts, f1, f2}};
       end
    end
end
end