function [x, ii] = solveAdjointEqsADI(eqs, eqs_p, adjVec, objk, system)


numVars = cellfun(@numval, eqs)';
cumVars = cumsum(numVars);
ii = [[1;cumVars(1:end-1)+1], cumVars];

if iscell(objk)
   objk = objk{:};
end
objk = cat(objk);

% Above CAT means '.jac' is a single element cell array.  Extract contents.
rhs  = -objk.jac{1}';
if ~isempty(adjVec)
    % If adjVec is not empty, we are not at the last timestep (first in the
    % adjoint recurrence formulation). This means that we subtract
    % the previous jacobian times the previous lagrange multiplier.
    % Previous here means at time t + 1.
    eqs_p = cat(eqs_p{:});

    % CAT means '.jac' is a single element cell array.
    rhs = rhs - eqs_p.jac{1}'*adjVec;
end
tic
if system.nonlinear.cpr
    [x, its, fl] = cprAdjoint(eqs, rhs, system, 'cprType', system.nonlinear.cprType, 'relTol', ...
                   system.nonlinear.cprRelTol, 'cprEllipticSolver', system.nonlinear.cprEllipticSolver);
else
    eqs = cat(eqs{:});

    % CAT means '.jac' is a single element cell array.
    x = eqs.jac{1}'\rhs;
end
tt = toc;
dispif(false, 'Lin. eq: %6.5f seconds, ', tt);
