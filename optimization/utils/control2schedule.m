function schedule = control2schedule(u, schedule, scaling)
% Convert control vector u to schedule 
nc = numel(schedule.control);
nw = numel(schedule.control(1).W);
[umin, umax] = deal(scaling.boxLims(:,1), scaling.boxLims(:,2));
c  = 0;
for cs = 1:nc
    for w = 1:nw
        c = c+1;
        schedule.control(cs).W(w).val = u(c) *(umax(w)-umin(w))+umin(w);
    end
end
end