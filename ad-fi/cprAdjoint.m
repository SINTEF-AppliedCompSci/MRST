function [lambda, its, fl] = cprAdjoint(eqs, rhs, system, varargin)
%% Description of the cpr preconditioner
% Let $A$ denote the full matrix given by the jacobian |cat(eqs{:})| and $d$ denote the
% right-hand side given by |rhs|.  We want to solve $A^t \lambda = d$.  From
% |eqsfiVO|, the columns of $A$  corresponding to the _main_ variables are $p$ (pressure), $s_w$
% (water saturation) and the last ones depends on the composition (oil, no gas or no oil,
% gas, or oil and gas). The CPR
% preconditioner of $A$ is applied after the following three preprocessing steps. 
%  
% # We switch the columns of $A$ so that the _main_ variables are $p$ (pressure), $s_w$
% (water saturation) and the last one corresponds to the remaining degree of freedom ($s_g$
% if oil and gas, $R_s$ if oil and no gas, $R_v$ if gas and no oil)
% # We remove directly the well variables.  
% # We assemble the pressure equations.
%
% Let $N_c$ be the number of cells and $N_w$ be the number of well variables. For each cell, we have
% $3$ equations for mass conservation (water, oil, gas). Thus, $A$ is a
% square matrix of dimension $N=3N_c+N_w$.
%
% The first preprocessing step (column switching) can be written as a right multiplication by a
% matrix which we denote $R$.
%
% The second preprocessing step is a block-wise Gaussian elimination, which can be written as a left
% multiplication by a matrix which we denote $L$.
%
% The third preprocessing step can also be written as a left multiplication by a matrix which we
% denote $L_p$.
%
% We denote by $\pi:R^N\to R^{3N_c}$ the projection mapping which selects
% from a _full_ vector of variables the vector with the _main_ variables. We denote by $\pi:R^N\to
% R^{N_c + N_w}$ the projection mapping which selects from the same _full_ vector the vector with
% the _secondary_ variables, the variables that are going to be eliminated in the second
% preprocessing step. By construction, we have
%
% $$\pi L A R \pi_c^t = 0.$$
%
%
% The preprocessing steps leads to the matrix $\tilde A$ of dimension $R^{3N_c\times 3N_c}$ given by
%
% $$ \tilde A = L_p \pi L A R\pi^t $$.
%
% The cpr preconditioner can be rewritten as a left multiplication by the matrix $M^{-1}$, which
% means that the iterative solver is applied to the matrix $M^{-1}\tilde A$.
%
%% The preprocessing equations for the adjoint solver.
% We want to solve
%
% $$A^t\lambda = d$$
%
% We multiply by $R^t$ and get
%
% $$R^tA^t\lambda = R^t d$$
%
% Let us define $\lambda_1$ as
%
% $$\lambda = L^t\lambda_1$$
%
% so that $\lambda_1$ is solution of
%
% $$R^tA^tL^t\lambda_1 = R^td.\quad\quad\textrm{(1)}$$
%
%
% We denote $\hat A=LAR$ and $\hat d=R^t d$ so that the last equation becomes
%
% $$\hat A^t\lambda_1 = \hat d.$$
%
% We decompose $\lambda_1$ into _main_ and _secondary_ variables:
%
% $$\lambda_1 = \pi^t\lambda_{1,1} + \pi_c^t\lambda_{1,2}$$
%
% We plug this into Equation (1) and get
%
% $$\hat A^t\pi^t\lambda_{1,1} + \hat A^t\pi_c^t\lambda_{1,2} = \hat d.$$
%
% We apply $\pi_c$ on both side and use the fact that $\pi_c\hat A \pi=0$ to obtain the equation
% for $\lambda_{1,2}$, i.e.,
%
% $$(\pi_c\hat{A}^t\pi_c^t)\,\lambda_{1,2} = \pi_c\hat{d} $$.
%
% The equation for $\lambda_{1,1}$ is then
%
% $$(\pi\hat{A}^t\pi^t)\,\lambda_{1,1} = \pi\hat{d} - \pi\hat{A}^t\pi_c^t\lambda_{1,2}.$$
%
% Let us introduce $\lambda_2$ as $\lambda_{1,1} = L_p^t \lambda_2$. Then the last equation
% yields
%
% $$\tilde A^t\,\lambda_2 = \pi\hat{d} - \pi\hat{A}^t\pi_c^t\lambda_{1,2}.$$
%
% To solve this equation for $\lambda_2$ (the right-hand side is at this point known), we use the
% same cpr preconditioner which is used for the forward equation (where we solve $Ax = b$).
%
%% Implementation
%

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

   opt = struct('cprType', 'colSum', 'cprEllipticSolver', @mldivide, 'relTol', 2e-2);
   opt = merge_options(opt, varargin{:});


   L = speye(numel(rhs));

   ne = numel(eqs);
   copt = system.cpr;
   active = max(copt.active);
   g = copt.gas;

   %%
   % *First preprocesing step.*
   %

   % If the equation were obtained via |eqsfiBlackOil|, we had to take
   % if ~isempty(g)
   %    % for unsaturated cells, switch respective columns of dx/ds and dx/drs
   %    unSat = logical(full(diag(eqs{g(2)}.jac{g(1)})));
   %    index = getEquationInxs(eqs);
   %    [eqs, rhs] = switchCols(eqs, rhs, index, g, unSat);
   % end
   %
   % However, the black oil equations are now by default assembled by ¦eqsfiVO| and in
   % this case, we do not have to reorder.

   %%
   % The variables |eqs| and |rhs| have now been update so that the full jacobi of eqs
   % corresponds to $AR$ and |rhs| to $R^td$.
   %

   full_eqs = cat(eqs{:});

   % CAT means '.jac' is a single element cell array.  Extract contents.
   fullA   = full_eqs.jac{1};

   %%
   % *Second preprocessing step*.
   %
   % We eliminate the well variables. We construct the matrices $L$, $\pi$ which is given by _proj_
   % and $\pi_c$ which is given by |proj_orth|.
   %

   n_elim = ne-active;

   proj = speye(numel(rhs));
   proj_orth = [];

   for i = 1:n_elim
      [eqs, L, proj, proj_orth] = elimVars(eqs, L, proj, proj_orth, active + 1);
   end

   fullA = L*fullA;

   ii = getEquationInxs(eqs);


   %%
   % *Third preprocessing step*.
   %
   % We assemble the pressure equations. We construct $L_p$.
   %

   [A, Lp, pInx] = getCPRSystemDiagonal(eqs, ii, opt);


   %%
   % We compute $\lambda_{1,2}$ which is given here by |lambda_wells|.
   %
   lambda_wells = (proj_orth*fullA'*proj_orth')\(proj_orth*rhs);

   %%
   % Now, we set |rhs| as $\pi\hat{d} - \pi\hat{A}^t\pi_c^t\lambda_{1,2}$:
   %
   rhs = proj*(rhs - fullA'*proj_orth'*lambda_wells);

   %%
   % We compute the preconditioner:
   %

   prec  = getTwoStagePrecAdj(A, pInx, opt);

   %%
   % We compute $\lambda_2$ which is here given by |lambda_cells|:
   %
   [lambda_cells, fl, relres, its] = gmres(A', rhs, [], opt.relTol, 20, prec);

   %%
   % |lambda|cells_ is updated to $\lambda_{1,1}$:
   %
   lambda_cells = Lp'*lambda_cells;

   %%
   % We compute $\lambda_1$ which is given here by |lambda|:
   %

   lambda = proj'*lambda_cells + proj_orth'*lambda_wells;

   %%
   % |lambda| is updated to $\lambda$, the solution of the adjoint equation:

   lambda = L'*lambda;

   if fl ~= 0
      warning('GMRES did not converge, Relative residual: %9.2e, error code: %2d', relres, fl);
   else
      fprintf('Used %d GMRES iterations\n', its(2));
   end
end

function eInx = getEquationInxs(eqs)
   numVars = cellfun(@numval, eqs)';
   cumVars = cumsum(numVars);
   eInx = [[1;cumVars(1:end-1)+1], cumVars];
end

%--------------------------------------------------------------------------
function [eqs, rhs]  = switchCols(eqs, rhs, ii, n, inx)
   index = find(inx); % inx is a logical;
   for k = 1:numel(eqs)
      % switch column for jacobian
      tmp = eqs{k}.jac{n(1)}(:, inx);
      eqs{k}.jac{n(1)}(:, inx) = eqs{k}.jac{n(2)}(:, inx);
      eqs{k}.jac{n(2)}(:, inx) = tmp;

   end
   % switch row for rhs
   tmp = rhs(ii(n(1), 1) + index - 1);
   rhs(ii(n(1), 1) + index - 1) = rhs(ii(n(2), 1) + index - 1);
   rhs(ii(n(2), 1) + index - 1) = tmp;
end

%--------------------------------------------------------------------------
function [eqs, L, proj, proj_orth] = elimVars(eqs, L, proj, proj_orth, n)
% remove eqs{n} by eliminating the corresponding unknown (i.e we invert eqs{n}.jac{n})

   solveInx = setdiff(1:numel(eqs), n);
   eqn      = eqs{n};

   index = getEquationInxs(eqs);
   L_iter = speye(index(end, end));

   inv_eqn_jac_n = inv(eqn.jac{n});
   for eqNum = solveInx
      for jacNum = solveInx
         eqs{eqNum}.jac{jacNum} = eqs{eqNum}.jac{jacNum} - eqs{eqNum}.jac{n}* ...
             (inv_eqn_jac_n*eqn.jac{jacNum});
      end
      L_iter(index(eqNum, 1) : index(eqNum, 2), index(n, 1):index(n, 2)) = L_iter(index(eqNum, 1) ...
                                                        : index(eqNum, 2), index(n, 1):index(n, ...
                                                        2)) - eqs{eqNum}.jac{n}*inv_eqn_jac_n;

      % Actually, we do not need to update the value of the residual for adjoint computation.
      eqs{eqNum}.val = eqs{eqNum}.val - eqs{eqNum}.jac{n}*(inv_eqn_jac_n*eqn.val);
   end

   eqs  = eqs(solveInx);
   for eqNum = 1:numel(eqs)
      eqs{eqNum}.jac = eqs{eqNum}.jac(solveInx);
   end


   L = proj'*L_iter*proj*L;
   if ~isempty(proj_orth)
      L = L +  proj_orth'*proj_orth;
   end

   proj_orth_step = speye(index(end, end));
   proj_orth_step = proj_orth_step(index(n,1):index(n,2), :);
   proj_orth = vertcat(proj_orth_step*proj, proj_orth);

   proj_step = speye(index(end, end));
   proj_step = proj_step(setdiff(1:index(end, end), index(n,1):index(n,2)), :);
   proj = proj_step*proj;


end


function [A, Lp, pInx] = getCPRSystemDiagonal(eqs, ii, opt)
   pInx = false(ii(end,end), 1);
   pInx(ii(1,1):ii(1,2)) = true;


   if strcmpi(opt.cprType, 'diag')
      cprFunc = @diag;
   elseif strcmpi(opt.cprType, 'colsum')
      cprFunc = @sum;
   end


   n = numel(eqs{1}.val);

   deqs = eqs;
   for k = 1:numel(eqs)
      for l = 1:numel(eqs{1}.jac)
         deqs{k}.jac{l} = spdiags(cprFunc(eqs{k}.jac{l})', 0, n,n);
      end
   end
   deqs = cat(deqs{:});

   % CAT means '.jac' is a single element cell array.  Extract contents.
   D = deqs.jac{1};

   eqs = cat(eqs{:});

   % CAT means '.jac' is a single element cell array.  Extract contents.
   A   = eqs.jac{1};

   Lp = inv(D);
   A = Lp*A;


end

function [A, b, pInx] = getCPRSystemBlockInvert(eqs, ii, opt, active)
   assert(active == 2 || active == 3);

   pInx = false(ii(end,end), 1);
   pInx(ii(1,1):ii(1,2)) = true;



   if strcmp(opt.cprType, 'diag')
      cprFunc = @diag;
   elseif strcmp(opt.cprType, 'colSum')
      cprFunc = @sum;
   end

   if active == 3
      % Solving for pressure and two other cell wise variables
      dss = cell(2,2);
      dps = cell(1,2);
      pI = 1;
      sI = [2,3];
      for ii = 1:2
         for jj = 1:2
            dss{ii,jj} = reshape(cprFunc(eqs{sI(ii)}.jac{sI(jj)}), [], 1);
         end
         dps{ii} = reshape(cprFunc(eqs{pI}.jac{sI(ii)}), [], 1);
      end

      n = numel(dss{1,1});
      inx = (1:n)';
      dtrmInv = 1./(dss{1,1}.*dss{2,2} - dss{2,1}.*dss{1,2});

      DssInv  = sparse([inx, inx  , inx+n, inx+n], ...
                       [inx, inx+n, inx  , inx+n], ...
                       [dss{2,2}, -dss{1,2}, -dss{2,1}, dss{1,1}].*(dtrmInv*[1 1 1 1]), ...
                       2*n, 2*n);


      Dps     = sparse([inx, inx], [inx, inx+n], [dps{1}, dps{2}], n, 2*n);
   else
      % A simpler two variable system
      pI = 1;
      sI = 2;
      dss = reshape(cprFunc(eqs{sI}.jac{sI}), [], 1);
      dps = reshape(cprFunc(eqs{pI}.jac{sI}), [], 1);

      n = numel(dss);
      inx = (1:n)';

      DssInv  = sparse(inx, ...
                       inx, ...
                       1./dss, ...
                       n,n);

      Dps     = sparse(inx, inx, dps, n, n);
   end


   M = Dps*DssInv;

   eqs = cat(eqs{:});

   % CAT means '.jac' is a single element cell array.  Extract contents.
   A   = eqs.jac{1};
   b   = -eqs.val;

   A(pInx, :) = A(pInx, :) - M*A(~pInx,:);
   b(pInx)    = b(pInx)    - M*b(~pInx);
end

%--------------------------------------------------------------------------

function prec = getTwoStagePrecAdj(A, pInx, opt)
   Ap     = A(pInx, pInx);

   setup.type = 'nofill';
   [L, U] = ilu(A', setup);

   prec   = @(r)applyTwoStagePrecAdj(r, A, L, U, Ap, pInx, opt);
end

%--------------------------------------------------------------------------

function x = applyTwoStagePrecAdj(r, A, L, U, Ap, pInx, opt)
   p = zeros(size(r));
   q = U\(L\r);
   s = r - A'*q;
   p(pInx) = opt.cprEllipticSolver(Ap', s(pInx));
   x = p + q;
end



