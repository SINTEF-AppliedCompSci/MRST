function pvec = getScaledParameterVector(setup, params)
values = applyFunction(@(p)p.getParameterValue(setup), params);
% scale values
u = cell(size(values));
for k = 1:numel(u)
    u{k} = params{k}.scale(values{k});
end
pvec = vertcat(u{:});  
end