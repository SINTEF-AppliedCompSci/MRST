function u = double2ADI(u, dummy)
assert(isa(u,'double'));
if isnumeric(dummy)
    % The dummy variable is also a double, we simply return without
    % changing it.
    return
end
nval  = numel(u);
jac  = cellfun(@(j) sparse(nval, size(j, 2)), ...
               dummy.jac, 'UniformOutput', false);
u = ADI(u,jac);
end
