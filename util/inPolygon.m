function inside = inPolygon(poly, ptn)
    tol = 1e-6;
    poly = poly - ptn;
    ind = [1:size(poly,1), 1];
    prev_side = 0;
    for e = 1:numel(ind) - 1
        v1 = poly(ind(e), :);
        v2 = poly(ind(e + 1), :);
        v1Xv2 = v2(1)*v1(2) - v1(1)*v2(2);
        current_side = sign(v1Xv2);
        if abs(v1Xv2) <= tol
            inside = false;
            return
        elseif prev_side == 0
            prev_side = current_side;
        elseif prev_side ~= current_side
            inside = false;
            return
        end
    end
    inside = true;
end