function [val, du, dx] = match_displacement_objfun(u, x, dirdofs, target, scaling)

   % make vector from displacements x
   assert(isvector(x)); % should have been made a vector before calling this function
   assert(numel(x) == numel(target));
   
   % discard values corresponding to Dirichlet boundaries.  These play no
   % role in the objective function evaluation here.
   x(dirdofs) = [];
   target(dirdofs) = [];
   
   [u, x] = initVariablesADI(u, x); %#ok

   error = x - target;

   val = -sum(error(:).^2)/numel(x.val) * scaling;
   
   du = val.jac{1}(:);
   dx = val.jac{2}(:);
   val = value(val);
   
end
