function linConstS   = setupConstraints(linConst, schedule, scaling)
% Setup linear constraints for scaled problem. Assumes linConst applies to
% all control steps.
[umin, umax] = deal(scaling.boxLims(:,1), scaling.boxLims(:,2));
D = diag(umax-umin);
A = linConst.A*D;
b = linConst.b - linConst.A*umin;
nc = numel(schedule.control);
Ar = repmat({A}, [1, nc]);
linConstS = struct('A', blkdiag(Ar{:}), ...
                   'b', repmat(b, [nc, 1]) );
end
