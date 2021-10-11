function [val, du, dx] = total_displacement_objfun(u, x)

   krull = load('u_010');
   
   SCALING = 1e8; % @@ how better to do this?
   
   % make vector from displacements x
   assert(isvector(x)); % should have been made a vector before calling this function

   [u, x] = initVariablesADI(u, x); %#ok

   offset = x - krull.dd;

   val = sum(offset(:).^2)/numel(x.val) * SCALING/100;
   
   %val = -sum(x(:).^2)/numel(x.val) * SCALING;
   
   du = val.jac{1}(:);
   dx = val.jac{2}(:);
   val = value(val);
   
end
