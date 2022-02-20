function D = SolveTOFEqsADI(eqs, state, W, computeTracer, linsolve)
%Solve the time of flight equations for solveStationaryPressure.
%
% SYNOPSIS:
%  D = SolveTOFEqsADI(eqs, state, W, true)
%
% DESCRIPTION:
%   Solves tof/tracer problems efficiently for solveStationaryPressure.
%   This function is primarily meant to be called from other functions and
%   the interface reflects this. For a more general time of flight solver,
%   see computeTOFandTracer.
%
% REQUIRED PARAMETERS:
%
%   eqs   - Residual equations for pressure, time of flight and tracer.
%
%   state - Reservoir state
%
%   W     - The wells for which tracers are to be computed.
%
%   computeTracer - Boolean indicating if tracers are required.
%
% RETURNS:
%   D     - Time of flight / tracer structure.

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
if nargin < 5
    linsolve = @mldivide;
end

[tof1, itracer, inj] = solveTOF(eqs, 'forward', state, W, computeTracer, linsolve);
[tof2, ptracer, prod] = solveTOF(eqs, 'backward', state, W, computeTracer, linsolve);

D.tof = [tof1, tof2];
D.itracer = itracer;
D.ptracer = ptracer;
D.inj = inj;
D.prod = prod;

if computeTracer
    [val, D.ppart] = max(D.ptracer,[],2);
    [val, D.ipart] = max(D.itracer,[],2); %#ok<*ASGLU>
else
    D.ppart = [];
    D.ipart = [];
end

end



function [tof, tracer, isInj] = solveTOF(eqs, direction, state, W, computeTracer, linsolve)

    if strcmpi(direction, 'forward')
        self = 4;
        other = 5;
        sgn = 1;
    else
        self = 5;
        other = 4;
        sgn = -1;
    end

    isInj = find(arrayfun(@(x) sgn*sum(x.flux) >= 0, state.wellSol)) .';

    e = eqs{self};
    rhs = e.val;
    %A = tofRobustFix(-e.jac{self});
    A = -e.jac{self};


    % Well flux and pressure influence on local system
    % adjust = e.jac{2}*sgn*vertcat(state.wellSol.flux) - e.jac{1}*state.pressure;

    % These are zero as long as initial guess is zero !!!
    adjust = e.jac{2}*vertcat(state.wellSol.flux) + e.jac{1}*state.pressure;


    % Well flux influence on the other time of flight system. We will use
    % the inverted source terms from the other system to drive tracer flow.
    adjust_tracer = eqs{other}.jac{2}*sgn*vertcat(state.wellSol.flux);

    rhs = rhs + adjust;
    if computeTracer
        rhs = [rhs, zeros(numel(rhs), numel(isInj))];
        for i = 1:numel(isInj)
            wc = W(isInj(i)).cells;
            tmp = 0*e.val;
            tmp(wc) = -sgn;

            tmp = tmp.*adjust_tracer;

            rhs(:, i+1) = tmp;
        end
    end
    sol = linsolve(A,rhs);
    sol(isnan(sol)) = deal(0);
    sol(isinf(sol)) = deal(max(sol(~isinf(sol))));

    tof = sol(:,1);
    tof(tof < 0) = 10*max(tof);
    if computeTracer
        tracer = max(sol(:, 2:end), 0);
    else
        tracer = [];
    end
end
