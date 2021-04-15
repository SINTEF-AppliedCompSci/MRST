function [pdu] = projectGradient(controls, du, varargin)
% Project gradient according to linear input constraints. Handles box-constraints and
% linear equality - and inequality-constraints based on iteratively applying the constraints
% to the gradient until convergence.
%
% SYNOPSIS:
%  [pGrad] = projectGradient(...
%    simRes, G, S, W, rock, fluid, schedule, controls, grad, objectiveFunction, varargin)
%
% DISCRIPTION:
%  Project gradient according to constraints iteratively until tollerance ConstTol is met
%  or max number of iterations MaxConstIts is met (returnes failure).
%  Constraints are applied in the following order:
%     i) Box const.
%    ii) Lin. InEq. const.
%   iii) Lin. Eq. const.
%
% PARAMETERS:
%   controls, grad
%
%   'pn'/pv - List of 'key'/value pairs defining optional parameters.  The
%             supported options are:
%          - ConstTol            : Relative contraint satisfaction tol
%          - MaxConstIts         : max number of its for constraint satisfaction
%          - VerboseLevel        : amount of output to screen
%
% RETURNS:
%   pGrad      - projected gradient
%
% SEE ALSO:

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


opt     = struct('ConstTol',    1e-3, ...
                 'MaxConstIts',  100, ...
                 'VerboseLevel',   0);
opt          = merge_options(opt, varargin{:});

%--------------------------------------------------------------------------
if opt.VerboseLevel >= 2, fprintf('\n********* Incorporating linear input constraints ***********\n'); end

numSteps = size(du, 2);
pdu      = du;
uInit    = [controls.well.values]';
normU    = norm(uInit, 'inf');

%--------------------------------------------------------------------------
tol = normU*opt.ConstTol/2;

minMax   =  vertcat(controls.well.minMax);
boxMin = minMax(:,1) + tol;
boxMax = minMax(:,2) - tol;

if ~isempty( controls.linEqConst )
    Aeq   = controls.linEqConst.A;   Beq   = controls.linEqConst.b;
else
    Aeq = []; Beq = [];
end
if ~isempty( controls.linIneqConst )
    Aineq = controls.linIneqConst.A; Bineq = controls.linIneqConst.b;
else
    Aineq = []; Bineq = [];
end

itCounts = zeros(numSteps, 1);

for step = 1 : numSteps
    uInitN = uInit(:, step);
    pGradN = du(:, step);
    numIneq = numel(Aineq);
    boxOK = false;
    if isempty(Aineq)
        ineqOK = true(numIneq, 1);
    else
        ineqOK = false(numIneq, 1);
    end
    if isempty(Aeq), eqOK = true; else eqOK = false; end

    itCount = 0;
    while (~boxOK || ~all(ineqOK) || ~eqOK) && (itCount < opt.MaxConstIts)
        itCount = itCount + 1;
        % 1. box constraints
        boxOK = false;
        pGradNNew = min( max( pGradN, boxMin-uInitN), boxMax-uInitN );
        if norm(pGradNNew-pGradN, 'inf') < tol
            boxOK = true;
        end
        pGradN = pGradNNew;
        % 2. linear ineq. const.
        for k = 1:numel(Aineq)
            a = Aineq{k};
            b = Bineq{k};
            c = (a' * (uInit + pGradN) - b)/(a'*a);
            if norm(a*c, 'inf') < tol % c either small or negative
                ineqOK(k) = true;
            else
                ineqOK(k) = false;
                pGradN = pGradN - a*c;
            end
        end
        % 3. linear eq. const.
        if ~isempty(Aeq)
            C = (Aeq*Aeq')\(Aeq*(uInitN + pGradN) - Beq);
            if abs( norm(Aeq'*C, 'inf') ) < tol
                eqOK = true;
            else
                eqOK = false;
                pGradN = pGradN - Aeq'*C;
            end
        end
    end
    itCounts(step) = itCount;
    pdu(:, step) = pGradN;
end
if any(itCounts == opt.MaxConstIts)
    warning('Some constraints may not be satisfied or feasible direction not found');
end
if opt.VerboseLevel >= 2, fprintf('Maximum number of its: %4d. ', max(itCounts)); end









