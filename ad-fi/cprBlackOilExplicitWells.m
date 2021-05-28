function [dx, its, fl] = cprBlackOilExplicitWells(eqs, varargin)
%Undocumented Utility Function

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

opt = struct('cprType'    ,'colSum', ...
             'ellipSolve' ,@mldivide, ...
             'relTol'     , 2e-2,...
             'pscale'     , 1);

opt = merge_options(opt, varargin{:});

% assume for now that eqs are ordered as {O, W, G, closure, wells}
%                     vars are odered as {pO, sW, sG, rs, BHP}

% for unsaturated cells, switch respective columns of dx/ds and dx/drs
unSat = logical(full(diag(eqs{4}.jac{3})));
eqs   = switchCols(eqs, [3 4], unSat);

% eliminate unknowns in position 4 (rs if saturated, sG otherwise)
[eqs, eq4] = elimVars(eqs, 4);

% Eliminate wells which also end up at position 4 once BHP and rs/sG has
% been eliminated
[eqs, eq_qWs] = elimVars(eqs, 4);
[eqs, eq_qOs] = elimVars(eqs, 4);
[eqs, eq_qGs] = elimVars(eqs, 4);

% eliminate BHP -unknowns (now at position 4)
[eqs, eqW] = elimVars(eqs, 4);

ii = getEquationInxs(eqs);

[A, b, pInx]   = getCPRSystem(eqs, ii, opt);

% Scale pressure variables
% pscale = 1e-7;



pscale = opt.pscale;
if pscale ~= 1;
    pind = ii(1,1):ii(1,2);
    A(:,pind) = A(:,pind) ./ pscale;
end

prec           = getTwoStagePrec(A, pInx, opt);


[dx, fl, relres, its] = gmres(A, b, [], opt.relTol, 20, prec);

if pscale ~= 1;
    dx(pind) = dx(pind) ./ pscale;
end

%[dx, fl, relres, its] = gmres(A, b, [], opt.relTol, 20);
dpO  = dx(ii(1,1):ii(1,2));
dsW  = dx(ii(2,1):ii(2,2));
dM1  = dx(ii(3,1):ii(3,2));  % sG/RS

% Recover eliminated variables. This is done in the reverse order of the
% elimination.
dBHP = recoverVars(eqW, 4,    {dpO, dsW, dM1});

dqGs = recoverVars(eq_qGs, 4, {dpO, dsW, dM1,                   dBHP});
dqOs = recoverVars(eq_qOs, 4, {dpO, dsW, dM1,             dqGs, dBHP});
dqWs = recoverVars(eq_qWs, 4, {dpO, dsW, dM1,       dqOs, dqGs, dBHP});

dM2  = recoverVars(eq4, 4,    {dpO, dsW, dM1, dqWs, dqOs, dqGs, dBHP});

%Assign dsG and dRS
dsG = ~unSat.*dM1  +  unSat.*dM2;
dRS =  unSat.*dM1  + ~unSat.*dM2;

dx = {dpO, dsW, dsG, dRS, dqWs, dqOs, dqGs, dBHP};
if fl ~= 0
    warning('GMRES did not converge, Relative residual: %9.2e, error code: %2d', relres, fl);
end
end

function eInx = getEquationInxs(eqs)
numVars = cellfun(@numval, eqs)';
cumVars = cumsum(numVars);
eInx = [[1;cumVars(1:end-1)+1], cumVars];
end

%--------------------------------------------------------------------------
function eqs   = switchCols(eqs, n, inx)
for k = 1:numel(eqs)
    tmp = eqs{k}.jac{n(1)}(:, inx);
    eqs{k}.jac{n(1)}(:, inx) = eqs{k}.jac{n(2)}(:, inx);
    eqs{k}.jac{n(2)}(:, inx) = tmp;
end
end

%--------------------------------------------------------------------------
function [eqs, eqn] = elimVars(eqs, n)
% eliminate set of unknowns nr n using equation n ()
solveInx = setdiff(1:numel(eqs), n);
eqn      = eqs{n};

for eqNum = solveInx
    for jacNum = solveInx
        eqs{eqNum}.jac{jacNum} = eqs{eqNum}.jac{jacNum} - eqs{eqNum}.jac{n}*(eqn.jac{n}\eqn.jac{jacNum});
    end
    eqs{eqNum}.val = eqs{eqNum}.val - eqs{eqNum}.jac{n}*(eqn.jac{n}\eqn.val);
end

eqs  = eqs(solveInx);
for eqNum = 1:numel(eqs)
    eqs{eqNum}.jac = eqs{eqNum}.jac(solveInx);
end

end
%--------------------------------------------------------------------------
function x = recoverVars(eq, n, sol)
% recover variables x at position n using solutions sol
solInx = [1:(n-1), (n+1):(numel(sol)+1)];
x = - eq.jac{n}\(eq.val);
for k  = 1:numel(solInx)
    x = x - eq.jac{n}\(eq.jac{solInx(k)}*sol{k});
end
end

function [A, b, pInx] = getCPRSystem(eqs, ii, opt)
% Assume vars are ordered pO, sW, sG/rs
% Eqs are ordered O, W, G
pInx = false(ii(end,end), 1);
pInx(ii(1,1):ii(1,2)) = true;

%sInx = ~pInx;

dss = cell(2,2);
dps = cell(1,2);
if strcmp(opt.cprType, 'diag')
    cprFunc = @diag;
elseif strcmp(opt.cprType, 'colSum')
    cprFunc = @sum;
end

pI = 1;
sI = [2,3];
for ii = 1:2
    for jj = 1:2
        dss{ii,jj} = reshape(cprFunc(eqs{sI(ii)}.jac{sI(jj)}), [], 1);
    end
    dps{ii} = reshape(cprFunc(eqs{pI}.jac{sI(ii)}), [], 1);
end

eqs = cat(eqs{:});

% CAT means '.jac' is a single element cell array.  Extract contents.
A   = eqs.jac{1};
b   = -eqs.val;

n = numel(dss{1,1});
inx = (1:n)';
dtrmInv = 1./(dss{1,1}.*dss{2,2} - dss{2,1}.*dss{1,2});
% fprintf('Cond det: %e \n', max(abs(dtrmInv))/min(abs(dtrmInv)));
% if any(dtrmInv<0)
%     fprintf(['Negative determinant in ', num2str(nnz(dtrmInv<0)), ' cells !!!\n' ])
% end
DssInv  = sparse([inx, inx  , inx+n, inx+n], ...
                 [inx, inx+n, inx  , inx+n], ...
                 [dss{2,2}, -dss{1,2}, -dss{2,1}, dss{1,1}].*(dtrmInv*[1 1 1 1]), ...
                 2*n, 2*n);


Dps     = sparse([inx, inx], [inx, inx+n], [dps{1}, dps{2}], n, 2*n);

M = Dps*DssInv;

A(pInx, :) = A(pInx, :) - M*A(~pInx,:);
b(pInx)    = b(pInx)    - M*b(~pInx);
end

%--------------------------------------------------------------------------

function prec = getTwoStagePrec(A, pInx, opt)
Ap     = A(pInx, pInx);
% setup.droptol = 1e-6;
%struct('type','ilutp','droptol',1e-5)
% 'nofill' for ilu(0)
setup.type = 'nofill';
[L, U] = ilu(A, setup);
% [L, U] = ilu(A);
prec   = @(r)applyTwoStagePrec(r, A, L, U, Ap, pInx, opt);
end

%--------------------------------------------------------------------------

function x = applyTwoStagePrec(r, A, L, U, Ap, pInx, opt)
x = zeros(size(r));
x(pInx) = opt.ellipSolve(Ap, r(pInx));
%x(pInx) = agmg(Ap, r(pInx));%, [], [], [], 1);
r = r - A*x;
x = x + U\(L\r);
end



