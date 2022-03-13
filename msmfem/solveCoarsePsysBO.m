function state = solveCoarsePsysBO(state, G, CG, p, rock, ...
                                   S, CS, fluid, p0, dt, varargin)
%Solve coarsened fine-scale well system (for Black Oil).
%
% SYNOPSIS:
%   [xr, xw] = solveCoarsePsysBO(xr, xw, G, CG, p, rock, S, CS, ...
%                                fluid, p0, dt)
%   [xr, xw] = solveCoarsePsysBO(xr, xw, G, CG, p, rock, S, CS, ...
%                                fluid, p0, dt, 'pn1', pv1, ...)
%
% REQUIRED PARAMETERS:
%   xr, xw   - Reservoir and well solution structures either properly
%              initialized from functions 'initResSol' and 'initWellSol'
%              respectively, or the results from a previous call to
%              function 'solveCoarsePsysBO'.
%
%   G, CG, p - Grid, coarse grid, and cell-to-block partition vector,
%              respectively.
%
%   rock     - Rock data structure.  Must contain valid field 'rock.perm'.
%
%   S, CS    - Linear system structure on fine grid (S) and coarse grid
%              (CS) as defined by functions 'computeMimeticIP' and
%              'generateCoarseSystem', respectively.
%
%   fluid    - Black Oil fluid object as defined by, e.g., function
%              'initBlackoilFluid'.
%
%   p0       - Vector, length G.cells.num, of cell pressures at previous
%              time step (not previous iteration of successive substitution
%              algorithm).
%
%   dt       - Time step size.
%
% OPTIONAL PARAMETERS:
%   wells    - Well structure as defined by functions 'addWell' and
%              'generateCoarseWellSystem'.  May be empty (i.e., W = [],
%              default value) which is interpreted as a model without any
%              wells.
%
%   bc       - Boundary condition structure as defined by function 'addBC'.
%              This structure accounts for all external boundary conditions
%              to the reservoir flow.  May be empty (i.e., bc = [], default
%              value) which is interpreted as all external no-flow
%              (homogeneous Neumann) conditions.
%
%   src      - Explicit source contributions as defined by function
%              'addSource'.  May be empty (i.e., src = [], default value)
%              which is interpreted as a reservoir model without explicit
%              sources.
%
%   LinSolve - Handle to linear system solver software to which the fully
%              assembled system of linear equations will be passed.
%              Assumed to support the syntax
%
%                        x = LinSolve(A, b)
%
%              in order to solve a system Ax=b of linear equations.
%              Default value: LinSolve = @mldivide (backslash).
%
% RETURNS:
%   xr, xw - Updated reservoir and well solution structures.  The reservoir
%            solution structure 'xr' will have added fields 'blockPressure'
%            and 'blockFlux' corresponding to the coarse grid pressure and
%            block fluxes, respectively.
%
% SEE ALSO:
%   `solveBlackOilWellSystem`.

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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


% This code implements one iteration of the coarsening scheme (12) from the
% paper
%
%    SPE 118993
%
%    A Multiscale Mixed Finite-Element Solver for
%    Three-Phase Black-Oil Flow
%
% by S. Krogstad, K.-A. Lie, H. M{\o}ll Nilsen, J.R. Natvig,
%    B. Skaflestad, and J.E. Aarnes, SINTEF
%
% As far as possible, the names of symbols in this code follow those of the
% cited reference.  In the main code, arrays with a suffix 'r' (e.g., 'Cr')
% refer to quantities in the reservoir while arrays with a suffix 'w'
% (e.g., 'Cw') refer to quantities in wells.  The sole exception is array
% 'wr' which denotes the index range (rows/columns) of matrix 'B'
% corresponding to wells.  Arrays with a prefix of 'd' (e.g., 'dF') refer
% to Dirichlet (pressure) conditions--either faces (dF) or the conditions
% values themselves (dC).  Also note that the matrices \Tilde{C} and
% \Tilde{C}_w of the paper are respectively termed 'Ctr' and 'Ctw' below.
%
% Subfunctions, e.g. 'nonlin_term_res', forego this convention as they
% typically focus on just one aspect of the global system (e.g., the
% reservoir part or the well part).
%
   opt = struct('bc', [], 'src', [], 'wells', [], ...
                'LinSolve', @mldivide,            ...
                'volume_correction', zeros([G.cells.num, 1]));
   opt = merge_options(opt, varargin{:});

   [nsub, sub] = get_subfaces(G, CG, CS);
   [I, J]      = coarsening_ops(G, CG, CS, sub, nsub);

   % Compute gravity terms.
   grav   = gravity();
   grav   = reshape(grav(1 : size(G.nodes.coords,2)), 1, []);

   nc     = G.cells.num;
   cf     = G.cells.faces(:,1);
   cellNo = rldecode(1 : nc, diff(G.cells.facePos), 2) .';
   sgn    = 2*double(G.faces.neighbors(cf, 1) == cellNo) - 1;

   % cfn  = cell-face normal = *outward* face normals on each cell.
   cfn    = bsxfun(@times, G.faces.normals(cf,:), sgn);

   % ng  == n' *     g for all cellfaces.
   % nKg == n' * K * g for all cellfaces.
   [K, row, col] = permTensor(rock, size(G.nodes.coords,2));

   assert (size(K,1) == G.cells.num, ...
          ['Permeability must be defined in active cells only.\n', ...
           'Got %d tensors, expected %d (== number of cells).'],   ...
           size(K,1), G.cells.num);

   ng     = cfn * reshape(grav, [], 1);
   nKg    = sum(cfn(:,row) .* bsxfun(@times, K(cellNo,:), grav(col)), 2);

   %-----------------------------------------------------------------------
   % Get updated pressure and saturation dependent quantities -------------
   %
   [c, rho, mu, u] = fluid.pvt(state.pressure, state.z);
   s     = bsxfun(@rdivide, u, u * ones([size(u,2), 1]));
   mob   = bsxfun(@rdivide, fluid.relperm(s), mu);

   % Parameters used in pressure equation.
   g2    = grav * grav.';      % Recall: 'grav' is a *ROW* vector.
   Lt    = sum(       mob                              , 2)      ;
   ct    = sum( s  .*  c                               , 2)      ;
   omega = sum(rho .* mob                              , 2) ./ Lt;
   zeta  = sum( c  .* mob                              , 2) ./ Lt;
   beta  = sum( c  .* mob .* bsxfun(@minus, omega, rho), 2)      ;

   Lti = sparse(1 : numel(cellNo), 1 : numel(cellNo), 1 ./ Lt(cellNo));

   [Bv, Phi]      = res_basis(G, CG, CS, p, Lt);
   [Psi_w, Phi_w] = well_basis(opt.wells, G, p, Lt, rho);
   flx_bas        = [Bv, -Psi_w];

   actB = CS.activeCellFaces;
   actF = CS.activeFaces;
   Psi  = S.BI * flx_bas;

   B  = (Lti * flx_bas).' * Psi;
   Vr = nonlin_term_res(faceFlux2cellFlux(G, state.flux), S.BI, Lt, ...
                        zeta, omega, beta, cellNo, ng, nKg);

   % Compute block matrices from:
   accum  = ct .* poreVolume(G, rock) ./ dt;
   P      = spdiags(-accum, 0, nc, nc);
   g_comp = p0 .* accum;

   g_grav = - G.cells.volumes .* beta .* omega .* g2;

   %  1) Wells
   wr = numel(actB) + 1 : size(B,1);  % Index range for wells in 'B'.
   if isfield(state, 'wellSol'),
      wellSol = state.wellSol;
   else
      wellSol = [];
   end
   [Bw, Cw, Dw, Ctw, fw, hw] = well_contrib(opt.wells, wellSol, ...
                                            G, p, Lt, omega, zeta, ...
                                            Psi(:,wr), Phi_w, P, Vr, I);

   if ~isempty(opt.wells),
      dFw = reshape(strcmpi({ opt.wells.type } .', 'bhp'), [], 1);
      dCw = reshape([ opt.wells(dFw).val ]               , [], 1);
   else
      dFw = logical([]);
      dCw = [];
   end

   %  2) Reservoir
   [Cr, Dr, Ctr, P] = res_contrib(Psi(:, 1:numel(actB)), ...
                                  Phi, CS.C(actB,:), CS.D(actB,actF), ...
                                  Vr, P, I);
   [fr, gr, hr, dFr, dCr] = sys_rhs(G, g_comp + g_grav + ...
                                    opt.volume_correction,    ...
                                    omega, Psi(:,1:numel(actB)), ...
                                    I, J, opt.src, opt.bc);

   dF = [dFr; dFw];
   dC = [dCr; dCw];

   % Add well contributions into mass matrix
   B(wr,wr) = B(wr,wr) + blkdiag(Bw{:});

   % Assemble final coarse system (block-matrix formulation).
   nc = size(Cr,2); nf = numel(actF) - sum(dF) + size(Dw,2);

   C  = [Cr ; Cw ];  D = blkdiag(Dr, Dw);  D = D(:,~dF);
   Ct = [Ctr; Ctw].';

   % System matrix.
   X = [  B ,                   C  ,            D  ; ...
         Ct ,                   P  , sparse(nc, nf); ...
         D.', sparse(size(D,2), nc +            nf)];

   ff = [fr; fw];
   gg =  gr     ;
   hh = [hr; hw];

   % Solve system.
   %
   if sprank(X) == size(X,1) - 1,
      % Force pressure zero level (for incompressible flows).
      szB = size(B);
      X(szB(1)+1, szB(1)+1) = 1;
   elseif sprank(X) < size(X,1),
      error('Singular coarse pressure matrix')
   end
   x = opt.LinSolve(X, [ff; gg; hh(~dF)]);

   lam      =  zeros(size(hh));
   flux     =  x(1 : numel(ff));
   pres     = -x(numel(ff) + 1 : numel(ff) + numel(gg));
   lam(~dF) =  x(numel(ff) + numel(gg) + 1 : end);
   lam( dF) = dC;

   % Assign output
   %  1) Reservoir
   state.blockPressure   = pres;
   state.blockFlux       = flux(1 : numel(actB));

   state.pressure(:)     = pres(p)                 + ...
                           Phi   * state.blockFlux + ...
                           Phi_w * (-flux(numel(actB) + 1 : end));
   state.flux(:)         = cellFlux2faceFlux(G, Psi * flux);

   if ~isfield(state, 'facePressure'),
      state.facePressure = zeros([G.faces.num, 1]);
   end

   state.facePressure(:) = accumarray(G.cells.faces(:,1), ...
                                      state.pressure(cellNo)) ./ ...
                           accumarray(G.cells.faces(:,1), 1);

   state.facePressure(sub) = rldecode(lam(1 : numel(actF)), nsub);

   %  2) Wells
   fluxW = flux(numel(actB) + 1 : end);
   lamW  = lam (numel(actF) + 1 : end);  assert (numel(lamW) == numel(opt.wells));

   nw = numel(opt.wells);
   i  = 0;
   for k = 1 : nw,
      rates = opt.wells(k).CS.rates;

      state.wellSol(k).flux     = - full(rates * fluxW(i + 1 : i + size(rates,2)));
      state.wellSol(k).pressure = lamW(k);

      i = i + size(rates,2);
   end
end


%--------------------------------------------------------------------------
% Helpers follow.
%--------------------------------------------------------------------------


function [nsub, sub] = get_subfaces(g, cg, cs)
% get_subfaces -- Extract fine-scale faces constituting all active coarse
%                 faces.
%
% SYNOPSIS:
%   [nsub, sub] = get_subfaces(G, CG, CS)
%
% PARAMETERS:
%   G, CG - Fine and coarse grid discretisation of the reservoir model,
%           respectively.
%
%   CS    - Coarse linear system structure as defined by function
%           'generateCoarseSystem'.
%
% RETURNS:
%   nsub, sub - The return values from function 'subFaces' compacted to
%               only account for active coarse faces (i.e. coarse faces
%               across which there is a non-zero flux basis function and,
%               therefore a possibility of non-zero flux).

   [nsub, sub] = subFaces(g, cg);
   sub_ix      = cumsum([0; nsub]);
   sub         = sub(mcolon(sub_ix(cs.activeFaces) + 1, ...
                            sub_ix(cs.activeFaces + 1)));
   nsub        = nsub(cs.activeFaces);
end

%--------------------------------------------------------------------------

function [I, J] = coarsening_ops(g, cg, cs, sub, nsub)
% coarsening_ops -- Compute coarsening operators from fine to coarse model.
%
% SYNOPSIS:
%   [I, J] = coarsening_ops(G, CG, CS, sub, nsub)
%
% PARAMETERS:
%   G, CG     - Fine and coarse grid discretisation of the reservoir model,
%               respectively.
%
%   CS        - Coarse linear system structure as defined by function
%               'generateCoarseSystem'.
%
%   sub, nsub - Packed representation of fine-scale faces constituting
%               (active) coarse faces.  Assumed to follow the conventions
%               of the return values from function 'subFaces'.
%
% RETURNS:
%  I - Coarse-to-fine cell scatter operator (such that I.' is the
%      fine-to-coarse block gather (sum) operator).
%
%  J - Coarse-to-fine face scatter operator for active coarse faces such
%      that J.' is the fine-to-coarse face gather (sum) operator.

   I = sparse(1:g.cells.num, cg.partition, 1);
   J = sparse(sub, rldecode(cs.activeFaces(:), nsub), ...
              1, g.faces.num, cg.faces.num);
   J = J(:, cs.activeFaces);
end

%--------------------------------------------------------------------------

function [C, D, Ct, P] = res_contrib(Psi, Phi, C, D, V, P, I)
% res_contrib -- Compute linear system contributions from reservoir.
%
% SYNPSIS:
%   [C, D, Ct, P] = res_contrib(Psi, Phi, C, D, V, P, I)
%
% PARAMETERS:
%   Psi, Phi   - Flux and pressure basis functions from reservoir blocks.
%
%   C, D, V, P - Fine-scale reservoir matrices C, D, V, and P,
%                resepectively.
%
%   I, J       - Coarsening operators as defined by function
%                'coarsening_ops'.
%
% RETURNS:
%   C, D, Ct, P - Coarse scale reservoir matrices C, D, \Tilde{C}, and P,
%                 respectively.

   assert (size(P,2) == size(P,1));

   PI = P * I;
   P  = I.' * PI;

   Ct = C  - (Psi.' * V) * I;
   Ct = Ct - Phi.' * PI;
end

%--------------------------------------------------------------------------

function [Bw, C, D, Ct, f, h] = well_contrib(W, xw, G, p, mob, omega, ...
                                             zeta, Psi, Phi, P, Vr, I)
% well_contrib -- Compute linear system contributions from wells.
%
% SYNOPSIS:
%   [B, C, D, Ct, f, h] = well_contrib(W, Psi, Phi, mob, C, P, Vr, Vw, I)
%
% PARAMETERS:
%   W        - Well structure as defined by functions 'addWell',
%              'assembleWellSystem', and 'generateCoarseWellSystem'.
%
%   Psi, Phi - Flux and pressure basis functions from well-to-block
%              connections.
%
%   mob      - Vector of (current) total mobility values.  One scalar value
%              for each cell in the underlying fine-scale model.
%
%   C, P     - Fine-scale cell-to-half-contact mapping (C) and accumulation
%              term (P).
%
%   Vr, Vw   - Fine-scale reservoir and well non-linear terms generated by
%              functions 'nonlin_term_res' and 'nonlin_term_well',
%              respectively.
%
%   I        - Coarse block-to-fine cell scatter operator (such that I.' is
%              the cell-to-coarse block gather operator).
%
% RETURNS:
%   B, C, D, Ct - Second term of matrix B_{22} and matrices C_w, D_w, and
%                 \Tilde{C}_w of equation (12) in the paper cited in the
%                 introduction.  Note that matrix B is stored in a cell
%                 array and must be passed to BLKDIAG (or similar) before
%                 being added into the global B_{22} matrix.
%
%   f, h        - Well contributions to linear system right hand side
%                 blocks 'f' and 'h' respectively.

   [Bw, C, D, f, h] = unpackWellSystemComponentsMS(W, G, p, mob, omega);

   C  = vertcat(C{:});
   D  = blkdiag(D{:});
   Ct = [];
   f  = vertcat(f{:});
   h  = vertcat(h{:});

   if ~isempty(W),
      nw = numel(W);
      Vw = cell([nw, 1]);

      for k = 1 : nw,
         wc    = W(k).cells;
         nc    = numel(wc);
         v     = zeta(wc) .* xw(k).flux ./ (mob(wc) .* W(k).WI);

         Vw{k} = W(k).CS.rates.' * sparse(1 : nc, wc, v, nc, G.cells.num);
      end

      Ct = (Psi.'*Vr + vertcat(Vw{:})) * I;
      Ct = C - Ct + (Phi.' * (P * I));
   end
end

%--------------------------------------------------------------------------

function [Bv, Phi] = res_basis(g, cg, cs, p, mob)
% res_basisP -- Assemble pressure basis functions for reservoir blocks.
%
% SYNOPSIS:
%   [Bv, Phi] = res_basis(G, CG, CS, p, mob)
%
% PARAMETERS:
%   G, CG - Fine-scale and coarse-scale discretisations (grids) of
%           reservoir model.
%
%   CS    - Coarse linear system structure as defined by function
%           'generateCoarseSystem'.
%
%   p     - Cell-to-coarse block partition mapping vector.
%
%   mob   - Vector of (current) total mobility values.  One scalar value
%           for each cell in the underlying fine-scale model.
%
% RETURNS:
%   Bv  - Reservoir flux basis function matrix (hybrid), B*\Psi.
%   Phi - Reservoir pressure basis functions represented as a sparse
%         matrix of size g.cells.num-by-nchc.

   [Bv, Phi] = basisMatrixHybrid(g, cg, cs);

   % Update pressure basis function matrix for effects of (modified?)
   % coarse block mobilities.  The average mobility is stored in the sixth
   % element of the input tuples (see evalBasisFunc for details).  We do
   % rely on the basis function average mobility to contain one or two
   % elements in the 'hybrid' case (two when the basis function is
   % supported in two distinct coarse blocks), and a single element in the
   % 'mixed' case.
   %
   % The update is a diagonal matrix containing the entries
   %    mob(block)_0 / mob(block)_now
   %
   mob0 = cellfun(@(x) x{6}, cs.basisP(cs.activeFaces), ...
                  'UniformOutput', false);
   mob  = accumarray(p, mob .* g.cells.volumes, [cg.cells.num, 1]) ./ ...
          accumarray(p,        g.cells.volumes, [cg.cells.num, 1]);

   blkno = rldecode(1 : cg.cells.num, diff(cg.cells.facePos), 2) .';
   mob   = vertcat(mob0{:}) ./ mob(blkno(cs.activeCellFaces));

   Phi = Phi * spdiags(mob, 0, size(Phi,2), size(Phi,2));
end

%--------------------------------------------------------------------------

function [Psi, Phi] = well_basis(w, g, p, mob, rho)
% well_basis -- Assemble flux and pressure basis functions for wells.
%
% SYNOPSIS:
%   [Psi, Phi] = well_basis(W, G, p, mob)
%
% PARAMETERS:
%   W, G - Fine-scale model Well and grid data structures respectively.
%
%   p    - Cell-to-coarse block partition mapping vector.
%
%   mob  - Vector of (current) total mobility values.  One scalar value for
%          each cell in the underlying fine-scale model.
%
% RETURNS:
%   Psi - Well flux basis functions represented as a single sparse matrix
%         of size nfhc-by-(sum_w (number of coarse blocks intersected by w)).
%
%   Phi - Well pressure basis functions represented as a single sparse
%         matrix of size g.cells.num-by-(sum_w (above)).

   t = cell([1, 5]);

   [t{:}, Psi, Phi] = unpackWellSystemComponentsMS(w, g, p, mob, rho);
end

%--------------------------------------------------------------------------

function V = nonlin_term_res(u, BI, mob, zeta, omega, ...
                             beta, cellNo, ng, nKg)
% nonlin_term_res -- Compute non-linear (velocity) term for reservoir.
%
% SYNOPSIS:
%   V = nonlin_term_res(u, G, S, fluid)
%
% PARAMETERS:
%   u     - Current outward cell fluxes (one scalar value for each
%           half-contact in underlying fine grid).
%
%   G, S  - Grid and linear system structures defining fine-scale model.
%
%   fluid - Black-oil fluid data structure as defined by functions
%           'readpvt' and 'pvt'.
%
% RETURNS:
%   V - Matrix, size nfhc-by-nfc, of values for the velocity term
%
%         V(u) = B * D_u * C * D_\zeta
%
%       though the values are computed a little differently due to the
%       non-availability of the matrix 'B'.
%
% NOTE:
%   This is a costly process due to the possibly (very) large dimension of
%   the matrix 'BI' (number of fine-scale half-contacts in the reservoir
%   model).

   n   = numel(cellNo);    i = 1 : n;

   BIV = bsxfun(@times, ng, beta(cellNo)) + ...
         zeta(cellNo) .* (u + omega(cellNo).*nKg);
   BIV = sparse(i, cellNo, BIV, n, cellNo(end));

   V   = (sparse(i, i, mob(cellNo), n, n) * BI) \ BIV;
end

%--------------------------------------------------------------------------

function [ff, gg, hh, dF, dC] = sys_rhs(g, gc, omega, Psi, I, J, src, bc)
% sys_rhs -- Evaluate coarse system right hand side contributions.
%
% SYNOPSIS:
%   [f, g, h, dF, dC] = sys_rhs(G, gc, omega, Psi, I, J, src, bc)
%
% PARAMETERS:
%   G     - Grid data structures.
%
%   gc    - Compressible RHS terms per fine scale cell.
%
%   omega - Accumulated phase densities \rho_i weighted by fractional flow
%           functions f_i -- i.e., omega = \sum_i \rho_i f_i.  One scalar
%           value for each cell in the discretised reservoir model, G.
%
%   Psi   - Flux basis functions from reservoir blocks.
%
%   I, J  - Coarsening operators as defined by function 'coarsening_ops'.
%
%   bc    - Boundary condition structure as defined by function 'addBC'.
%           This structure accounts for all external boundary conditions to
%           the reservoir flow.  May be empty (i.e., bc = struct([])) which
%           is interpreted as all external no-flow (homogeneous Neumann)
%           conditions.
%
%   src   - Explicit source contributions as defined by function
%           'addSource'.  May be empty (i.e., src = struct([])) which is
%           interpreted as a reservoir model without explicit sources.
%
% RETURNS:
%   f, g, h - Coarse-scale direct (block) right hand side contributions as
%             expected by the linear system solvers such as
%             'schurComplement'.
%
%   dF, dC  - Coarse-scale packed Dirichlet/pressure condition structure.
%             The faces in 'dF' and the values in 'dC'.  May be used to
%             eliminate known face pressures from the linear system before
%             calling a system solver (e.g., 'schurComplement').
%
% SEE ALSO:
%   `computePressureRHS`.

   % Evaluate fine-scale system right hand side contributions.
   %
   [ff, gg, hh, grav, dF, dC] = computePressureRHS(g, omega, bc, src);
   gg = gg + gc;

   % Accumulate fine-scale boundary conditions to coarse scale.
   %
   ff = Psi.' * (ff + grav);
   gg =  I .' * gg;
   hh  = J .' * hh;

   % Account for coarse scale Dirichlet boundary conditions.
   %
   dC2 = zeros(size(dF));  dC2(dF) = dC;
   dF  = J.' * dF;
   dC  = J.' * dC2 ./ dF;
   dF  = logical(dF);
   dC  = dC(dF);
end
