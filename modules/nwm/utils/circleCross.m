function p = circleCross(x1, y1, r1, x2, y2, r2)
% Compute intersction points of two circiles
    if y1 ~= y2
        d  = sqrt( (x1-x2)^2+(y1-y2)^2 );
        k2 = -(x1-x2)/(y1-y2);
        b2 = (r1^2-r2^2-x1^2+x2^2-y1^2+y2^2)/(2*(y2-y1));

        if abs(r2-r1) < d && d < r2+r1
            delta = -b2^2 + r2^2 + k2^2 *r2^2 - 2 *b2* k2* x2 - ...
                k2^2* x2^2 + 2*b2*y2 + 2*k2*x2*y2 - y2^2;
            xx1 = (-b2* k2 + x2 + k2 *y2 - sqrt(delta))/(1 + k2^2);
            yy1 = k2*xx1 + b2;
            xx2 = (-b2* k2 + x2 + k2 *y2 + sqrt(delta))/(1 + k2^2);
            yy2 = k2*xx2 + b2;
            p = [xx1,yy1; xx2,yy2];
        else
            p = nan(2,2);
        end
    else
        p0 = [x1, y1; x2, y2];
        R = [r1; r2];
        [~, idx] = sort(p0(:,1));
        p0 = p0(idx, :);
        R = R(idx);
        [x1, x2, y1, y2, r1, r2] = deal(p0(1,1), p0(2,1), p0(1,2), p0(2,2),...
            R(1), R(2));
        d = abs(x2-x1);
        delta = (-d+r1-r2)*(-d-r1+r2)*(-d+r1+r2)*(d+r1+r2);
        if delta > 0
            a = sqrt(delta) / d;
            dx = (d^2 - r2^2 + r1^2)/(2*d);
            p = [x1 + dx, y1+a/2; x1 + dx, y1-a/2];
        else
            p = nan(2,2);
        end
    end
end