function [mult, VW, report] = getPolymerShearMultiplier(model, ...
    VW0, muWmult)
% Compute the flux multiplier due to polymer shear thinning/thickening
%
% SYNOPSIS:
%   mult           = getPolymerShearMultiplier(model, VW0, muWmult)
%   [mult,VW]      = getPolymerShearMultiplier(...)
%   [mult,VW,iter] = getPolymerShearMultiplier(...)
%
% DESCRIPTION:
%   The viscoisty of a polymer solution may be changed by the shear rate,
%   which is related to the flow velocity. This function returns a flux
%   multiplier mult, such that
%     vWsh = vW.*mult,   and,   vPsh = vP.*mult,
%   where vWsh and vPsh are the shear modified fluxes for pure water and
%   water with polymer, respectively.
%
%   The input water face velocity is computed as
%     VW = vW/(poro*area)
%   where poro is the average porosity between neighboring cells, and area
%   is the face area.
%
% REQUIRED PARAMETERS:
%   model    - Model structure
%   VW0      - Water velocity on faces without shear thinning
%   muWmult  - Viscosity multiplier on faces
%   iter     - Number of non-linear iterations
%
% RETURNS:
%   mult     - Flux multiplier (reciprocal of the viscosity multiplier)
%   VW       - Water velocity as solution of nonlinear shear equation

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

t = tic;
f = model.fluid;

% Compute the flow velocity
VW0     = abs(double(VW0));
muWmult = double(muWmult);

% Setup equation for the shear modified velocity. This is an equation on
% the faces of the domain.
eqn = @(VW) VW .* (1 + f.plyshearMult(VW).*(muWmult-1)) - VW0.*muWmult;

% Solve for shear modified velocity
[VW, report]  = solveNonlinearEqn(eqn, VW0);

% Compute velocity multiplier from solution
% Note that this is the reciprocal of the viscosity multiplier
mult = muWmult ./ ( 1 + f.plyshearMult(VW).*(muWmult-1) );

% Add timing to report
report.Time = toc(t);

end



%--------------------------------------------------------------------------
% HELPER FUNCTIONS
%--------------------------------------------------------------------------

function [x, report] = solveNonlinearEqn(fun, x0)
% Simple Newton solver to find the root of the nonlienar shear equation

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
    % Tolerence not met. Issue warning.
    warning('Polymer shear convergence failure. Residual=%1.2e', ...
        maxit, resnorm);
end

x = double(x);

report.Iterations = iter;
report.Residual   = resnorm;
report.Converged  = resnorm <= abstol;
end
