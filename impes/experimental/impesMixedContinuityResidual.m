function [F, Fp, Fv] = impesMixedContinuityResidual(pv, dt, fluid, mob, ...
                                                    v, press, mass, W, OP)
%Form (black-oil) continuity equation contributions to (IMPES) mixed system
%
% SYNOPSIS:
%   [F, Fp, Fv] = impesMixedContinuityResidual(pv, dt, fluid, mob, ...
%                                              v, press, mass, OP)
%
% DESCRIPTION:
%   This function forms contributions to a linear system of equations of
%   mixed type that may subsequently be used either directly to recover
%   pressures and fluxes or, at the caller's discretion, in an algebraic
%   reduction of MsMFE type to a smaller system of equations in preparation
%   of a linear system solve.
%
%   Specifically, this function forms the blocks 'Fv', 'Fp', and 'F' of the
%   (block) system of linear equations.
%
%       [  B    C  ] [ v ]   [ f ]
%       [ Fv   -Fp ] [-p ] = [-F ]
%
%   which defines discrete pressure and flux values in a mixed
%   discretisation.  It is the caller's responsibility to form the matrix
%   blocks 'B', 'C', and 'f' before attempting a linear system solve.  The
%   function is primarily intended in a multiscale solution strategy.
%
%   This function uses a two-point flux approximation (TPFA) method with
%   minimal memory consumption within the constraints of operating on a
%   fully unstructured polyhedral grid structure.  Moreover, we attempt to
%   minimise the number of fluid matrix evaluations.
%
% PARAMETERS:
%   pv    - Pore volume.  One positive scalar for each grid cell.  Use
%           function 'poreVolume' to compute this value.
%
%   dt    - Time step size (IMPES/pressure step).
%
%   fluid - Compressible fluid object.
%
%   mob   - Phase-mobilities defined per connection/interface.  Function
%           'tpfaUpwindStateVars' may be used to compute this value based
%           on an upwind strategy.
%
%   v     - Total flux.  One scalar value for each interface.  Specifically,
%           NUMEL(v)==SIZE(mob,1) is assumed to hold.
%
%   press - Data structure from which to extract pressure values.  Assumed
%           to maintain at least the following two fields:
%              - cell  --  Scalar pressure values per cell
%              - face  --  Scalar pressure values per interface
%
%   mass  - Data structure from which to extract component masses or,
%           rather, the 'z' variable (component volume per pore volume).
%           Assumed to maintain at least the following two fields:
%              - cell  --  Component volume values per cell
%              - face  --  Component volume values per interface
%
%   OP    - TPFA connection operators as defined by function
%           'tpfaConnectionOperators'.
%
% NOTE:
%   The function 'tpfaUpwindStateVars' constitutes a possible definition of
%   the 'mass.face' state variable.
%
% RETURNS:
%   F  - IMPES residual at the phase-space point (v,press,mass).  Vector of
%        size nc-by-1 with 'nc' being NUMEL(press.cell).
%
%   Fp - Derivative of 'F' with respect to (cell) pressure at the point
%        (v,press,mass).  Diagonal, sparse matrix of size nc-by-nc.
%
%   Fv - Derivative of 'F' with respect to total velocity ('v') at the
%        point (v,press,mass).  Sparse matrix of size nc-by-NUMEL(v).
%
% SEE ALSO:
%   impesTPFA, tpfaConnectionOperators, tpfaUpwindStateVars, poreVolume.

%{
Copyright 2009-2015 SINTEF ICT, Applied Mathematics.

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


   [Ac, dAc] = fluid_matrix(fluid, press.cell, mass.cell);
   Af        = fluid_matrix(fluid, press.face, mass.face);

   nc = numel(pv);
   np = size(mass.cell, 2);

   n  = nc * np;
   if size(Ac, 1) > n,
      i   = 1 : n;
      Ac  = Ac (i, i);
      dAc = dAc(i, i);
   end

   luAc = factorise(Ac);

   % DISREGARD GRAVITY FOR THE TIME BEING
   Alam = mmultiply(Af, bsxfun(@rdivide, mob, sum(mob, 2)));

   cellno = [OP.cellno ; OP.wcellno];  connno = [OP.connno ; OP.wconnno];
   sgn    = [OP.sgn    ; OP.wsgn   ];

   z = bsxfun(@times, pv, mass.cell(1:nc, :)) - ...
       sparse(cellno, 1 : numel(cellno), dt .* sgn .* v(connno)) ...
       *                                                         ...
       Alam(connno, :);

   mat_bal = solve_single_sys(luAc, z);

   Fc  = pv - sum(mat_bal, 2);

   Fpc = mmultiply(dAc, mat_bal);
   Fpc = sparse(1 : numel(pv), 1 : numel(pv), ...
                sum(solve_single_sys(luAc, Fpc), 2));

   Fvc = solve_multiple_sys(luAc, Alam(connno,:), OP);
   Fvc = sparse(cellno, connno, dt .* sgn .* sum(Fvc, 2));

   if ~isempty(OP.wsgn),
      nperf_tot = sum(cellfun('prodofsize', { W.cells }));

      Ap   = Af(end - nperf_tot*np + 1 : end, ...
                end - nperf_tot*np + 1 : end);
      pmob = mob(end - nperf_tot + 1 : end, :);
      pflx = v(OP.wconnno);
      nf   = size(mob, 1) - nperf_tot;

      [Fw, Fpw, Fvw] = well_equation(W, nf, pflx, Ap, pmob);
   else
      [Fw, Fpw, Fvw] = deal([]);
   end

   F  =        [Fc  ; dt .* Fw ];
   Fp = blkdiag(Fpc , dt .* Fpw);
   Fv =        [Fvc ; dt .* Fvw];
end

%--------------------------------------------------------------------------

function [A, dA, varargout] = fluid_matrix(fluid, p, z)
   %c    rho  mu  u  R  B  A  dA
   [rho, rho, A,  A, A, A, A, dA] = fluid.pvt(p, z);  %#ok

   if nargout > 2,
      varargout{1} = rho;
   end
end

%--------------------------------------------------------------------------

function [F, Fp, Fv] = well_equation(W, nf, v, Ap, pmob)
   nperf   = reshape(cellfun('prodofsize', { W.cells }), [], 1);
   is_bhp  = reshape(strcmpi('bhp'       , { W.type  }), [], 1);
   is_resv = reshape(strcmpi('resv'      , { W.type  }), [], 1);

   if ~all(is_bhp),
      pickw = ~is_bhp;
      nw    = sum(pickw);

      pickp = rldecode(pickw, nperf);
      resvp = rldecode(is_resv, nperf);

      fac = mmultiply(Ap, bsxfun(@rdivide, pmob, sum(pmob, 2)));
      fac = sum(fac .* rldecode(vertcat(W.compi), nperf), 2);
      fac(resvp) = 1;

      fac = fac(pickp);

      wno = rldecode(cumsum(double(pickw)), nperf);
      F   = accumarray(wno(pickp), fac .* v(pickp)) - [ W(pickw).val ].';
      Fp  = sparse(nw, nw);
      Fv  = sparse(wno(pickp), nf + find(pickp), fac, ...
                   nw        , nf + sum (nperf)     );
   else
      [F, Fp, Fv] = deal([]);
   end
end

%--------------------------------------------------------------------------

function v = mmultiply(A, x)
   np = size(x, 2);
   v  = reshape(A * reshape(x .', [], 1), np, []) .';
end

%--------------------------------------------------------------------------

function luA = factorise(A)
   assert (issparse(A), 'Internal error in ''%s''.', mfilename);

   [luA.L, luA.U, luA.P, luA.Q, luA.R] = lu(A);
end

%--------------------------------------------------------------------------

function x = solve_single_sys(luA, b)
   np = size(b, 2);
   x  = solve_sys(luA, reshape(b .', [], 1));
   x  = reshape(x, np, []) .';
end

%--------------------------------------------------------------------------

function x = solve_multiple_sys(luA, b, OP)
   np = size(b, 2);

   B = sparse(OP.I, OP.J, reshape(b .', [], 1));

   X = solve_sys(luA, B);

   x = reshape(full(X(sub2ind(size(X), OP.I, OP.J))), np, []) .';
end

%--------------------------------------------------------------------------

function x = solve_sys(luA, b)
   x = luA.Q * (luA.U \ (luA.L \ (luA.P * (luA.R \ b))));
end
