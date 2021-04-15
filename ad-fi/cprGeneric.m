function [dx, its, fl] = cprGeneric(eqs, system, varargin)
% A generic CPR preconditioner for implicit systems.
% Currently assumes that
% - The first variable is the pressure
% - That any cellwise variables to be solved for are next.
% - That any variables which can be eliminated before solving comes after
% the cellwise variables.
% - Well closure equations are last.

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
                'pscale'     , system.pscale, ...
                'eqScale'    , [1, 1, 1]);

   opt = merge_options(opt, varargin{:});
   ne = numel(eqs);

   copt = system.cpr;
   %scale eqs:
   for eqn = 1:numel(copt.active)
       eqs{eqn} = eqs{eqn}.*opt.eqScale(eqn);
   end

   active = max(copt.active);
   g = copt.gas;
   if and(~isempty(g), numel(g)==2)
      % for unsaturated cells, switch respective columns of dx/ds and dx/drs
      unSat = logical(full(diag(eqs{g(2)}.jac{g(1)})));
      eqs   = switchCols(eqs, g, unSat);
   end

   n_elim = ne-active;
   eliminated = cell(n_elim,1);

   for i = 1:n_elim
      [eqs, eq] = elimVars(eqs, active + 1);
      eliminated{i} = eq;
   end

   ii = getEquationInxs(eqs);

   if system.nonlinear.cprBlockInvert
      [A, b, pInx] = getCPRSystemBlockInvert(eqs, ii, opt, active);
   else
      [A, b, pInx] = getCPRSystemDRS(eqs, ii, opt, active);
      %[A, b, pInx] = getCPRSystemDiagonal(eqs, ii, opt);
   end

   % Scale pressure variables
   pscale = opt.pscale;
   if pscale ~= 1;
      pind = ii(1,1):ii(1,2);
      A(:,pind) = A(:,pind) ./ pscale;
   end

   prec = getTwoStagePrec(A, pInx, opt);

   [dX, fl, relres, its] = gmres(A, b, [], opt.relTol, 40, prec);
   if pscale ~= 1;
      dX(pind) = dX(pind) ./ pscale;
   end

   dx = cell(ne, 1);
   for i = 1:active
      dx{i} = dX(ii(i,1):ii(i,2));
   end

   % Recover eliminated variables. This is done in the reverse order of the
   % elimination.

   for i = n_elim:-1:1
      dVal = recoverVars(eliminated{i}, active + 1, {dx{[1:active (active + i + 1 : ne)]}}); %#ok<CCAT1>
      dx{active + i} = dVal;
   end


   %Assign dsG and dRS
   if and(~isempty(g), numel(g)==2)
      dM1 = dx{g(1)};
      dM2 = dx{g(2)};
      dx{g(1)} = ~unSat.*dM1  +  unSat.*dM2;
      dx{g(2)} =  unSat.*dM1  + ~unSat.*dM2;
   end

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
         if(~(prod(size(eqs{eqNum}.jac{jacNum}))==0) && ~( prod(size(eqn.jac{jacNum}))==0) )
            eqs{eqNum}.jac{jacNum} = eqs{eqNum}.jac{jacNum} - eqs{eqNum}.jac{n}*(eqn.jac{n}\eqn.jac{jacNum});
         end
      end
      if(~isempty(eqs{eqNum}.val) && ~isempty(eqn.val) )
         eqs{eqNum}.val = eqs{eqNum}.val - eqs{eqNum}.jac{n}*(eqn.jac{n}\eqn.val);
      end

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
      if(~isempty(x) && ~isempty(sol{k}))
         x = x - eq.jac{n}\(eq.jac{solInx(k)}*sol{k});
      end
   end
end

function [A, b, pInx] = getCPRSystemDRS(eqs, ii, opt, active)
if active == 2
    pInx = false(ii(end,end), 1);
    pInx(ii(1,1):ii(1,2)) = true;

    edd = .1; % 0<= edd <= 1
    epsilon = 0;
    n = numel(eqs{1}.val);
    l1 = zeros(n,2);
    for k = 1:2
        dj  = diag(eqs{k}.jac{1});
        sod = ( sum(abs(eqs{k}.jac{1}))-abs(dj)' )';
        l1(:,k) = dj./sod > edd;
    end
    % check for zero elems of d
    l21 = zeros(n,1); l22 = ones(n,1);
    ix = find(l1(:,1)==0);
    if ~isempty(ix)
        chk = l1(ix,2)>0;
        l1(ix(~chk),1) = 1;
        l21(ix(chk)) = 1;
        l22(ix(chk)) = 0;
    end
    i1 = (1:n)'; i2 = i1+n;

    L = sparse([i1, i1, i2 , i2 ], ...
               [i1, i2, i1 , i2 ], ...
               [l1,     l21, l22], ...
               2*n, 2*n);

    eqs = cat(eqs{:});
    A   = eqs.jac{1};
    b   = -eqs.val;

    A = L*A;
    b = L*b;
elseif active == 3
    pInx = false(ii(end,end), 1);
    pInx(ii(1,1):ii(1,2)) = true;

    edd = .5; % 0<= edd <= 1
    epsilon = 0;
    n = numel(eqs{1}.val);
    l1 = zeros(n,3);
    for k = 1:3
        dj  = diag(eqs{k}.jac{1});
        sod = ( sum(abs(eqs{k}.jac{1}))-abs(dj)' )';
        l1(:,k) = dj./sod > edd;
    end
    % check for zero elems of d
    l21 = zeros(n,1); l22 = ones(n,1);
    l31 = zeros(n,1); l33 = ones(n,1);
    ix = find(l1(:,1)==0);
    if ~isempty(ix)
        l12x = l1(ix,2);
        l13x = l1(ix,3);
        bothZero = (l12x+l13x)==0;
        l1(ix(bothZero),1) = 1;
        g23 = l12x>=l13x;
        l21(ix( and(~bothZero, g23) ),1) = 1;
        l22(ix( and(~bothZero, g23) ),1) = 0;
        l31(ix( and(~bothZero, ~g23) ),1) = 1;
        l33(ix( and(~bothZero, ~g23) ),1) = 0;
    end

    % %remove 'weak couplings'
    % tt = false(n,2);
    % dj  =  abs(diag(eqs{1}.jac{1}));
    % for k = 1:2
    %     sod = sum(abs(eqs{1}.jac{k}),2);
    %     tt(:,k) = (sod < epsilon*dj);
    % end
    % ix  = find(tt(:,1));
    % chk = and( or(l1(ix,1),l1(ix,3)), l22(ix) );
    % l1(ix(chk),2) = 0;
    %
    % ix  = find(tt(:,2));
    % chk = and( or(l1(ix,1),l1(ix,2)), l33(ix) );
    % l1(ix(chk),3) = 0;




    i1 = (1:n)'; i2 = i1+n; i3 = i2+n;

    L = sparse([i1, i1, i1, i2 , i2 , i3 , i3 ], ...
        [i1, i2, i3, i1 , i2 , i1 , i3 ], ...
        [l1,         l21, l22, l31, l33], ...
        3*n, 3*n);

    eqs = cat(eqs{:});
    A   = eqs.jac{1};
    b   = -eqs.val;

    A = L*A;
    b = L*b;
    % disp([sum(l1), sum(tt)])
end
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

function prec = getTwoStagePrec(A, pInx, opt)
   Ap     = A(pInx, pInx);

   setup.type = 'nofill';
   [L, U] = ilu(A, setup);

   prec   = @(r)applyTwoStagePrec(r, A, L, U, Ap, pInx, opt);
end

%--------------------------------------------------------------------------

function x = applyTwoStagePrec(r, A, L, U, Ap, pInx, opt)
   x = zeros(size(r));
   x(pInx) = opt.ellipSolve(Ap, r(pInx));

   r = r - A*x;
   x = x + U\(L\r);
end
