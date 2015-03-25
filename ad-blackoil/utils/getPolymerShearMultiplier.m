function mult = getPolymerShearMultiplier(model, vW0, muWmultf)
% Compute the velocity multiplier caused by polymer shear thinning
% 
% The viscoisty of a polymer solution may be changed by the shear rate,
% which is related to the flow velocity. This function returns a multiplier
% M, such that vWsh = vW.*M, and vPsh = vP.*M, where vWsh and vPsh are the
% shear modified velocities of pure water and water with polymer,
% respectively.
% 

f = model.fluid;

% Compute the flow velocity
vW0      = abs(double(vW0));
muWmultf = double(muWmultf);

% Setup equation for the shear modified velocity. This is an equation on
% the faces of the domain.
eqn = @(v) v .* (1 + f.plyshearMult(v).*(muWmultf-1)) - vW0.*muWmultf;

% Solve for shear modified velocity
v   = solveNonlinearEqn(eqn, vW0);

% Compute multiplier from solution
mult = muWmultf ./ ( 1 + f.plyshearMult(v).*(muWmultf-1) );

end



%--------------------------------------------------------------------------
% HELPER FUNCTIONS
%--------------------------------------------------------------------------

function x = solveNonlinearEqn(fun, x0)

% Hardcoded settings
maxit   = 30;
abstol  = 1.e-15;

% Initialize
iter = 0;
x    = x0;

% Newton loop
while iter <= maxit
   
   x   = initVariablesADI(x);
   eqs = fun(x);
   
   resnorm = norm(double(eqs), 'inf');
   if resnorm < abstol
      break
   end
   
   J   = -eqs.jac{1};
   dx  = J\eqs.val;
   x   = x.val + dx;
   
   iter = iter + 1;
   
end

if resnorm > abstol
  error('Polymer shear convergence failure. Residual=%1.2e', ...
     maxit, resnorm);
end

x = double(x);

end



