function mult = getPolymerShearMultiplier(model, VW0, muWmultf)
% Compute the velocity multiplier due to polymer shear thinning/thickening
% 
% SYNOPSIS:
%   mult = getPolymerShearMultiplier(model, vW0, muWmultf)
%
% PARAMETERS:
%   model    - Model structure
%   vW0      - Water velocity on faces without shear thinning (see below)
%   muWmultf - Viscosity multiplier on faces
%
% The viscoisty of a polymer solution may be changed by the shear rate,
% which is related to the flow velocity. This function returns a velocity 
% multiplier mult, such that
%   vWsh = vW.*mult,   and,   vPsh = vP.*mult,
% where vWsh and vPsh are the shear modified velocities of pure water and 
% water with polymer, respectively.
% 
% The water given water face velocity is computed as
%   VW = vW/(poro*area)
% where poro is the average porosity between neighboring cells, and area is
% the face area.
% 

f = model.fluid;

% Compute the flow velocity
VW0      = abs(double(VW0));
muWmultf = double(muWmultf);

% Setup equation for the shear modified velocity. This is an equation on
% the faces of the domain.
eqn = @(VW) VW .* (1 + f.plyshearMult(VW).*(muWmultf-1)) - VW0.*muWmultf;

% Solve for shear modified velocity
VW  = solveNonlinearEqn(eqn, VW0);

% Compute velocity multiplier from solution
% Note that this is the reciprocal of the viscosity multiplier
mult = muWmultf ./ ( 1 + f.plyshearMult(VW).*(muWmultf-1) );

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



