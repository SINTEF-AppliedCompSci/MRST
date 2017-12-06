function normal = normalFromPoints(p)
    assert(size(p,1)>=3)
    p1 = p(1, :);
    p2 = p(2,:);
    p_mean = mean(p, 1);
    v1 = p2 - p1;
    v2 = p_mean - p1;
    
    normal = cross(v1, v2);
    if sum(abs(normal))<1e-8
        normal = normalFromPoints(p(:, 2:end));
        return
    end
    normal = normal / sqrt(sum(normal.^2));
end