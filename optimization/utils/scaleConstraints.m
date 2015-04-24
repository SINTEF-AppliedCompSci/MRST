function linConstS   = scaleConstraints(linConst, scaling)
[umin, umax] = deal(scaling.boxLims(:,1), scaling.boxLims(:,2));
D = diag(umax-umin);
linConstS = struct('A'     , linConst.A*D, ...
                   'b'     , linConst.b - linConst.A*umin, ...
                   'scaled', true);
end
