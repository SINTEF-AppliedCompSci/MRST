function u = double2ADI(u, dummy)
assert(isa(u,'double'));
nval  = numel(u);
jac  = cellfun(@(j) sparse(nval, size(j, 2)), ...
               dummy.jac, 'UniformOutput', false);
u = ADI(u,jac);
end
