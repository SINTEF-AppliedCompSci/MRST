function u = schedule2control(schedule, scaling)
nc = numel(schedule.control);
[umin, umax] = deal(scaling.boxLims(:,1), scaling.boxLims(:,2));
u = cell(nc, 1);
for c = 1:nc
    ui   = vertcat(schedule.control(c).W(:).val);
    u{c} = (ui-umin)./(umax-umin);
end
u = vertcat(u{:});
end

