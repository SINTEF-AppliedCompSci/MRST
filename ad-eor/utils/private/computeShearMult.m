function v = computeShearMult(fluid, Vw, muWMultf)
% Compute the shear multiplier by solving EQ 52.12 in the ECLIPSE Technical
% Description.

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

% Vw should be the absolute value?
    f = fluid;
    % The solution of the shear multipler will be performed through an
    % iterative non-linear solution of the EQ. 52.12 in TD.
    % V ( 1+(P-1)M(V) ) / P = Vw;
    % P is the muWmultf, which is from PLYVISC

    % give the initial guess of the Vsh
    Vsh = Vw;

    Vsh = initVariablesADI(Vsh);

    plyshearMult = f.plyshearMult;

    shFunc = @(x) x.*(1+(muWMultf-1.).*plyshearMult(x))-muWMultf.*Vw;
    eqs = shFunc(Vsh);

    resnorm = norm(value(eqs), 'inf');
    iter = 0;
    maxit = 30;
    abstol = 1.e-15;

    while (resnorm > abstol) && (iter <= maxit)

      J = eqs.jac{1};
      dVsh = -(J \ eqs.val);
      Vsh.val = Vsh.val + dVsh;

      eqs = shFunc(Vsh);
      resnorm = norm(value(eqs), 'inf');

      iter = iter + 1;

    end

    if (iter >= maxit) && (resnorm > abstol)
        error('Convergence failure within %d iterations\nFinal residual = %.8e', maxit, resnorm);
    else
        M = plyshearMult(Vsh.val);
        v = (1 + (muWMultf - 1.).* M) ./ muWMultf;
    end
end
