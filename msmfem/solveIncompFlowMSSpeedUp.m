function state = solveIncompFlowMSSpeedUp(state, G, CG, p, S, CS, fluid, varargin)
%Solve (multiscale) pressure system assuming that some values are precomputed.
%
% SYNOPSIS:
%   state = solveIncompFlowMS(state, G, CG, p, S, CS, fluid)
%   state = solveIncompFlowMS(state, G, CG, p, S, CS, fluid, ...
%                             'pn1', pv1, ...)
%
% DESCRIPTION:
%   This function assembles and solves a (block) system of linear equations
%   defining interface fluxes and cell and interface pressures at the next
%   time step in a sequential splitting scheme for the reservoir simulation
%   problem defined by Darcy's law and a given set of external influences
%   (wells, sources, and boundary conditions).
%
%   NB:
%   The difference between this function and the normal 'solveIncompFlowMS'
%   is that we demand that the CS structure contains precomputed values
%   from the basis functions in fields 'Bv' and 'Phi' while CG contains
%   precomupted values for the subFaces in fields 'nsub' and 'sub'.  This
%   fields are added to CS and CG by the function 'speedUpMS'
%
% REQUIRED PARAMETERS:
%   state  - Reservoir and well solution structure either properly
%            initialized from functions 'initResSol' and 'initWellSol'
%            respectively, or the results from a previous call to function
%            'solveIncompFlowMS' and, possibly, a transport solver such as
%            function 'implicitTransport'.
%
%   G, CG, p -
%            Grid, coarse grid, and cell-to-block partition vector,
%            respectively.
%
%   S, CS  - Linear system structure on fine grid (S) and coarse grid
%            (CS) as defined by functions 'computeMimeticIP' and
%            'generateCoarseSystem', respectively.
%
%   fluid  - Fluid object as defined by function 'initSimpleFluid'.
%
% OPTIONAL PARAMETERS:
%   wells  - well linear system structure (W) as defined by functions
%            'addWell', 'assembleWellSystem', and
%            'assembleCoarseWellSystem'.  An empty well linear system
%            structure (i.e., W=struct([])) is supported.  The resulting
%            model is interpreted as not containing any well driving
%            forces.
%
%   bc     - Boundary condition structure as defined by function 'addBC'.
%            This structure accounts for all external boundary conditions
%            to the reservoir flow.  May be empty (i.e., bc = struct([]))
%            which is interpreted as all external no-flow (homogeneous
%            Neumann) conditions.
%
%   src    - Explicit source contributions as defined by function
%            'addSource'.  May be empty (i.e., src = struct([])) which is
%            interpreted as a reservoir model without explicit sources.
%
%   rhs    - Supply system right-hand side 'b' directly.  Overrides
%            internally constructed system right-hand side.  Must be a
%            three-element cell array, the elements of which are correctly
%            sized for the block system component to be replaced.
%
%            NOTE: This is a special purpose option for use by code which
%            needs to modify the system of linear equations directly, e.g.,
%            the 'adjoint' code.
%
%   Solver - Which solver mode function 'solveIncompFlowMS' should employ
%            in assembling and solving the block system of linear
%            equations.
%            String.  Default value: Solver = 'hybrid'.
%
%            Supported values are:
%              - 'hybrid' --
%                 Assemble and solve hybrid system for interface pressures.
%                 System is eventually solved by Schur complement reduction
%                 and back substitution.  The hybrid solver can not be
%                 employed if any of the basis functions in the coarse
%                 system were generated with a positive overlap.
%
%               - 'mixed' --
%                  Assemble and solve a mixed system for cell pressures and
%                  interface fluxes.
%
%   LinSolve -
%            Handle to linear system solver software to which the fully
%            assembled system of linear equations will be passed.  Assumed
%            to support the syntax
%
%                        x = LinSolve(A, b)
%
%            in order to solve a system Ax=b of linear equations.
%            Default value: LinSolve = @mldivide (backslash).
%
%   MatrixOutput -
%            Whether or not to return the final system matrix 'A' to the
%            caller of function 'solveIncompFlowMS'.
%            Logical.  Default value: MatrixOutput = FALSE.
%
% RETURNS:
%   xr - Reservoir solution structure with new values for the fields:
%          - cellPressure -- Pressure values for all cells in the
%                            discretised reservoir model, 'G'.
%          - facePressure -- Pressure values for all interfaces in the
%                            discretised reservoir model, 'G'.
%          - cellFlux     -- Outgoing flux across each local interface for
%                            all cells in the model.
%          - faceFlux     -- Flux across global interfaces corresponding to
%                            the rows of 'G.faces.neighbors'.
%          - A            -- System matrix.  Only returned if specifically
%                            requested by setting option 'MatrixOutput'.
%
%   xw - Well solution structure array, one element for each well in the
%        model, with new values for the fields:
%           - flux     -- Perforation fluxes through all perforations for
%                         corresponding well.  The fluxes are interpreted
%                         as injection fluxes, meaning positive values
%                         correspond to injection into reservoir while
%                         negative values mean production/extraction out of
%                         reservoir.
%           - pressure -- Well pressure.
%
% NOTE:
%   If there are no external influences, i.e., if all of the structures
%   'W', 'bc', and 'src' are empty and there are no effects of gravity,
%   then the input values 'xr' and 'xw' are returned unchanged and a
%   warning is printed in the command window.  This warning is printed with
%   message ID
%
%           'solveIncompFlowMS:DrivingForce:Missing'
%
% SEE ALSO:
%   `generateCoarseSystem`, `computeMimeticIP`, `addBC`, `addSource`, `addWell`,
%   `initSimpleFluid`, `initResSol`, `initWellSol`, `solveIncompFlow`.

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


   opt = struct('bc', [], 'src', [], 'wells', [], 'rhs', [], ...
                'Solver',       'hybrid',                    ...
                'LinSolve',     @mldivide,                   ...
                'MatrixOutput', false);
   opt = merge_options(opt, varargin{:});

   do_solve = check_input(G, CG, CS, state, opt);

   if do_solve,
      [A, b, dF, dC, Phi, Psi] = build_system(state, G, CG, p, S, CS, ...
                                              fluid, opt);
      solver = pick_solver(CG, CS, dF, dC, opt);
      x      = solver(A, b);
      state  = pack_solution(state, G, CG, p, CS, x{1:3}, Phi, Psi, opt);

      if opt.MatrixOutput, state.A = x{4}; end
   end
end

%--------------------------------------------------------------------------
% Helpers follow.
%--------------------------------------------------------------------------

function do_solve = check_input(G, CG, CS, state, opt)
   if isfield(state, 'wellSol'),
      xw = state.wellSol;
   else
      xw = [];
   end
   if ~isempty(opt.wells) && isempty(xw)
      error(id('wellSol:Empty'), ...
           ['Well solution structure must be initialized ', ...
            'by call to function initWellSol.\n']);
   end

   pressure_bc = ~isempty(opt.bc) && ...
                 any(strcmpi('pressure', opt.bc.type));
   well_bhp    = ~isempty(opt.wells) && ...
                 any(strcmpi('bhp', { opt.wells.type }));
   if ~(pressure_bc || well_bhp),
      sum_rates = 0;
      if ~isempty(opt.wells), sum_rates = sum(vertcat(opt.wells.val)); end
      if ~isempty(opt.bc), sum_rates = sum_rates + sum(opt.bc.value);  end
      if (abs(sum_rates) > eps)
         error(id('MassBalance:NotFulfilled'), ...
              ['Well rates and flux bc must sum up to 0 \n', ...
               'when there are no bhp constrained wells or pressure bc.\n']);
      end
   end

   if ~isempty(opt.bc),
      checkBoundaryBasis(G, CG, CS, opt.bc);
   end

   g_vec = gravity();
   grav  = norm(g_vec(1 : size(G.nodes.coords,2))) > 0;

   % We assemble and solve a system if there are any external forces or
   % sources (e.g., gravity, boundary conditions, explicit source terms,
   % injection/production wells or, in the case of the adjoint method, if
   % the caller supplied the right hand side directly).
   do_solve = grav || ~all([isempty(opt.bc),    ...
                            isempty(opt.src),   ...
                            isempty(opt.wells), ...
                            isempty(opt.rhs)]);
   if ~do_solve,
      warning(id('DrivingForce:Missing'),                      ...
             ['No external driving forces present in model--', ...
              'state remains unchanged.\n']);
   end

   % Assemble (and solve) system even in absence of external driving forces
   % if the caller requested 'MatrixOutput'.
   do_solve = do_solve || opt.MatrixOutput;
end

%--------------------------------------------------------------------------

function [A, b, dF, dC, Phi, Psi] = ...
         build_system(state, G, CG, p, S, CS, fluid, opt)
   [nsub, sub] = get_subfaces(G, CG, CS);
   [I, J]      = coarsening_ops(G, CG, p, CS, sub, nsub);

   [totmob, omega] = dynamic_quantities(state, fluid);

   % Build system components B, C, D, f, g, and h for final hybrid system
   %
   %    [  B   C   D  ] [  v ]     [ f ]
   %    [  C'  0   0  ] [ -p ]  =  [ g ] .
   %    [  D'  0   0  ] [ cp ]     [ h ]
   %
   % The components are represented as one-dimensional cell arrays as
   %
   %    A = { B, C, D }
   %    b = { f, g, h }
   %
   % both of which are assembled separately from reservoir contributions
   % (syscomp_res) and well contributions (syscomp_well).  We assemble the
   % global components by concatenating the well components onto the
   % reservoir components.  This assembly is a result of treating wells as
   % faces/contacts.
   %
   % Specifically, we assemble the components as
   %
   %    A{1} = B =   Br   +    Bw    % A{1}(wr,wr) = ...
   %    A{2} = C = [ Cr   ;    Cw ]  % VERTCAT
   %    A{3} = D = [ Dr, 0; 0, Dw ]  % BLKDIAG
   %
   %    b{1} = f = [ fr   ;    fw ]  % VERTCAT
   %    b{2} = g = [ gr   ;       ]  % No direct contr. to block from well
   %    b{3} = h = [ hr   ;    hw ]  % VERTCAT
   %
   [Ar,               Psi_r, Phi_r] = syscomp_res(G, CG, CS, p, totmob, opt);
   [Aw, bw, dFw, dCw, Psi_w, Phi_w] = syscomp_wells(opt.wells, G, p, ...
                                                    totmob, omega);

   Psi    = S.BI * [Psi_r, -Psi_w];
   Phi    =        [Phi_r,  Phi_w];

   [br{1:3}, dFr, dCr] = sys_rhs(G, omega, Psi(:,1:size(Psi_r,2)), ...
                                 I, J, opt.src, opt.bc);

   cellno = rldecode(1 : G.cells.num, diff(G.cells.facePos), 2) .';
   ncf    = numel(cellno);
   Ar{1}  = [Psi_r, -Psi_w].' * ...
            sparse(1 : ncf, 1 : ncf, 1 ./ totmob(cellno)) * Psi;

   A    = cell([1, 3]);              b    = cell([1, 3]);
   A{1} = Ar{1}                   ;  b{1} = vertcat(br{1}, bw{1}{:});
   A{2} = vertcat(Ar{2}, Aw{2}{:});  b{2} =         br{2}           ;
   A{3} = blkdiag(Ar{3}, Aw{3}{:});  b{3} = vertcat(br{3}, bw{2}{:});

   wr          = size(Psi_r,2) + (1 : size(Psi_w,2));
   A{1}(wr,wr) = A{1}(wr,wr) + blkdiag(Aw{1}{:});
   if strcmpi(opt.Solver, 'hybrid'),
      A{1} = invert(A{1}, A{2});
   end

   % If rhs is supplied by user, replace b constructed above by given rhs
   if ~isempty(opt.rhs),
      assert (iscell(opt.rhs));
      assert (numel(b) == numel(opt.rhs));
      assert (all(cellfun(@(u,v) all(size(u) == size(v)), b, opt.rhs)));
      b = opt.rhs;
   end

   % Finally, assemble a list of prescribed (contact) pressures, i.e., both
   % known face pressure and known bottom hole pressures.  These will be
   % eliminated from the linear system prior to calling one of the system
   % solvers, and entered directly into the solution component 'cp'.
   %
   dF   = [dFr; dFw];
   dC   = [dCr; dCw];
end

%--------------------------------------------------------------------------

function [A, Bv, Phi] = syscomp_res(g, cg, cs, p, mob, opt)
   % Form hybrid or mixed basis function matrices depending on user
   % requests and basis function suitability (cannot employ hybrid
   % formulation if any basis function were generated with positive
   % overlap).
   %

   if ~isfield(cs, {'Bv', 'Phi'})
      if strcmpi(opt.Solver, 'hybrid'),
         if any(cellfun(@(x) x{7} > 0, cs.basisP(cs.activeFaces))),
            error('solveIncompFlowMS:Solver:NotApplicable', '%s\n%s',    ...
               'The ''hybrid'' coarse solver is not applicable when', ...
               '    any basis function is generated on an extended domain.');
         end

         [Bv, Phi] = basisMatrixHybrid(g, cg, cs);


      else
         [Bv, Phi] = basisMatrixMixed (g, cg, cs);
      end

   else  % SPEED-UP HACK
      Bv = cs.Bv;
      Phi = cs.Phi;
   end

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
   Phi   = Phi * spdiags(mob, 0, size(Phi,2), size(Phi,2));

   % Extract the 'C' and 'D' matrices for active coarse faces.
   A    = cell([1, 3]);
   A{2} = cs.C(cs.activeCellFaces,        :      );
   A{3} = cs.D(cs.activeCellFaces, cs.activeFaces);
end

%--------------------------------------------------------------------------

function [A, b, dF, dC, Psi, Phi] = syscomp_wells(w, g, p, mob, omega)
   A = cell([1, 3]);
   b = cell([1, 2]);

   [A{:}, b{:}, Psi, Phi] = unpackWellSystemComponentsMS(w, g, p, ...
                                                         mob, omega);

   if isempty(w),
      dF = logical([]);
      dC =         [] ;
   else
      dF = reshape(strcmp({ w.type } .', 'bhp'), [], 1);
      dC = reshape([ w(dF).val ]               , [], 1);
   end
end

%--------------------------------------------------------------------------

function [nsub, sub] = get_subfaces(G, CG, CS)
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
if ~isfield(CG, {'nsub', 'sub'});
   [nsub, sub] = subFaces(G, CG);
   sub_ix      = cumsum([0; nsub]);
   sub         = sub(mcolon(sub_ix(CS.activeFaces) + 1, ...
                            sub_ix(CS.activeFaces + 1)));
   nsub        = nsub(CS.activeFaces);

else  % SPEED-UP HACK
   nsub = CG.nsub;
   sub  = CG.sub;
end
end

%--------------------------------------------------------------------------

function [I, J] = coarsening_ops(G, CG, p, CS, sub, nsub)
% coarsening_ops -- Compute coarsening operators from fine to coarse model.
%
% SYNOPSIS:
%   [I, J] = coarsening_ops(G, CG, p, CS, sub, nsub)
%
% PARAMETERS:
%   G, CG, p  - Fine and coarse grid discretisation of the reservoir model,
%               as well as cell-to-coarse block partition vector.
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

   I = sparse(1 : G.cells.num, p, 1, G.cells.num, CG.cells.num);
   J = sparse(sub, rldecode(CS.activeFaces(:), nsub), ...
              1, G.faces.num, CG.faces.num);
   J = J(:, CS.activeFaces);
end

%--------------------------------------------------------------------------

function [ff, gg, hh, dF, dC] = sys_rhs(G, omega, Psi, I, J, src, bc)
% sys_rhs -- Evaluate coarse system right hand side contributions.
%
% SYNOPSIS:
%   [f, G, h, dF, dC] = sys_rhs(G, omega, Psi, I, J, src, bc)
%
% PARAMETERS:
%   G     - Grid data structure.
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
   [ff, gg, hh, grav, dF, dC] = computePressureRHS(G, omega, bc, src);

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

%--------------------------------------------------------------------------

function BI = invert(B, C)
   % Pre-allocate storage for sparse input arrays.
   %
   n = sum(sum(C > 0) .^ 2);
   i = zeros([n, 1]);
   j = zeros([n, 1]);
   s = zeros([n, 1]);

   % Compute inverse mass matrix values by inverting local matrices
   % (one dense matrix for each coarse block).
   %
   ix = 0;
   for k = 1 : size(C,2),
      r   = reshape(find(C(:,k)), [], 1);
      nF  = numel(r);
      nF2 = nF^2;

      i(ix + 1 : ix + nF2) = repmat(r  , [nF, 1]);
      j(ix + 1 : ix + nF2) = repmat(r.', [nF, 1]);
      s(ix + 1 : ix + nF2) = inv( B(r,r) );

      ix = ix + nF2;
   end

   % Assemble final inverse mass matrix.
   %
   BI = sparse(i, j, s, size(B,1), size(B,2));
end

%--------------------------------------------------------------------------

function s = id(s)
   s = ['solveIncompFlowMS:', s];
end

%--------------------------------------------------------------------------

function solver = pick_solver(CG, CS, dF, dC, opt)
   regul = ~any(dF);  % Set zero level if no prescribed pressure values.

   %-----------------------------------------------------------------------
   % Enter prescribed pressures into final solution vector ----------------
   %
   lam     = zeros(size(dF));
   lam(dF) = dC;

   switch lower(opt.Solver),
      case 'hybrid', solver = @(A,b) solve_hybrid(A, b, lam, dF, regul, ...
                                                  opt.LinSolve);
      case 'mixed' , solver = @(A,b) solve_mixed(A, b, lam, CG, CS, dF, ...
                                                 regul, opt);
      otherwise,
         error(id('SolverMode:NotSupported'), ...
               'Solver mode ''%s'' is not supported.', opt.solver);
   end
end

%--------------------------------------------------------------------------

function x = solve_hybrid(A, b, lam, dF, regul, LinSolve)
   A{3} = A{3}(:,~dF);
   b{3} = b{3}(  ~dF);
   [x{1:4}] = schurComplementSymm(A{:}, b{:},          ...
                                  'Regularize', regul, ...
                                  'LinSolve', LinSolve);
   lam(~dF) = x{3};
   x{3}     = lam;
end

%--------------------------------------------------------------------------

function x = solve_mixed(A, b, lam, CG, CS, dF, regul, opt)
   nF = neumann_faces   (CG, dF, CS.activeFaces, opt);
   Do = oriented_mapping(CG, CS, opt);

   % Reduce size of system by only taking into account the 'Neumann'
   % faces.
   %
   A{3} = A{3}(:,nF);
   b{3} = b{3}(  nF);

   % We have already formed the 'mixed' B matrix in function
   % 'build_system'.  Thus, pass option 'MixedB' to linear solver.
   %
   [x{1:4}] = mixedSymm(A{:}, b{:}, Do, 'Regularize', regul, ...
                        'MixedB', true, 'LinSolve', opt.LinSolve);

   % Put known contact pressures into solution vector.  Note: These
   % contact pressures are only defined non-trivially on the external
   % boundary of the domain.
   %
   lam(nF) = x{3};
   x{3}    = lam;
end

%--------------------------------------------------------------------------

function state = pack_solution(state, G, CG, p, CS, flux, pres, ...
                               lam, Phi, Psi, opt)
   actB = CS.activeCellFaces;
   actF = CS.activeFaces;

   %-----------------------------------------------------------------------
   % Package reservoir (coarse) solution structure ------------------------
   %
   switch lower(opt.Solver),
      case 'hybrid',
         nhf      = numel(actB);
         flux_sgn = rldecode([1, -1], [nhf, size(Phi,2) - nhf], 2) .';

         state.blockFlux = flux(1 : numel(actB));
         fluxW           = flux(numel(actB) + 1 : end);
      case 'mixed',
         nf       = numel(actF);
         flux_sgn = rldecode([1, -1], [nf, size(Phi,2) - nf], 2) .';

         orient       = get_orientation(CG, CS);
         renum        = zeros([CG.faces.num, 1]);
         renum(actF)  = 1 : numel(actF);

         state.blockFlux = orient .* flux(renum(CG.cells.faces(actB,1)));
         fluxW           = flux(numel(actF) + 1 : end);
   end

   state.blockPressure = pres;
   state.pressure(:)   = pres(p) + Phi * (flux .* flux_sgn);

   if strcmpi(opt.Solver, 'hybrid'),
      % The 'mixed' solver does not generate reasonable coarse interface
      % pressure on internal coarse faces.  Therefore, only insert coarse
      % interface pressures computed by the solver in the 'hybrid' case.
      %
      if ~isfield(state, 'facePressure'),
         state.facePressure = zeros([G.faces.num, 1]);
      end
      state.facePressure(:) = nan;

      [nsub, sub] = get_subfaces(G, CG, CS);
      state.facePressure(sub) = rldecode(lam(1 : numel(actF)), nsub);
   end

   state.flux(:)     = cellFlux2faceFlux(G, Psi * flux);

   %-----------------------------------------------------------------------
   % Recover fine scale well fluxes and pressures. ------------------------
   %
   w    = opt.wells;
   lamW = lam(numel(actF) + 1 : end);  assert (numel(lamW) == numel(w));

   nw = numel(w);
   i  = 0;
   for k = 1 : nw,
      rates = w(k).CS.rates;

      state.wellSol(k).flux(:)  = - full(rates * fluxW(i + 1 : i + size(rates,2)));
      state.wellSol(k).pressure = lamW(k);

      i = i + size(rates,2);
   end
end

%--------------------------------------------------------------------------

function nF = neumann_faces(CG, dF, actF, opt)
% Determine the 'Neumann faces' (nF) (i.e., the faces for which
% face/contact pressures will be determined in the 'mixed' and 'TPFA'
% cases) of the model (g) according to the following rules:
%
%   - No internal faces are counted amongst the 'Neumann faces'.
%   - An external face is a 'Neumann face' unless there is a prescribed
%     pressure value (i.e., a Dirichlet condition) associated to the face.
%   - Wells controlled by rate-constraints are treated as Neumann contacts.
%

   nF = false(numel(actF),1); % No int.
   % All external active faces ... :
   nF(any(CG.faces.neighbors(actF,:) == 0, 2)) = true;
   nF(dF(1:numel(actF))) = false; % ... except Dirichlet cond.

   if ~isempty(opt.wells),
      % Additionally include all rate-constrained wells.
      nF = [nF; strcmp('rate', { opt.wells.type } .')];
   end
end

%--------------------------------------------------------------------------

function orient = get_orientation(cg, cs)
   blkno  = rldecode(1 : cg.cells.num, diff(cg.cells.facePos), 2) .';
   c1     = cg.faces.neighbors(cg.cells.faces(:,1), 1);
   orient = 2*double((blkno == c1)) - 1;
   orient = orient(cs.activeCellFaces);
end

%--------------------------------------------------------------------------

function Do = oriented_mapping(CG, CS, opt)
% 'Do' maps face fluxes to half-face fluxes.  This matrix is used to form
% the reduced, mixed system of linear equtions in linear system solver
% function 'mixedSymm'.

   orient = get_orientation(CG, CS);

   if ~isempty(opt.wells),
      ws = [ opt.wells.CS ];
      dw = diag(vertcat(ws.D));
      ow = ones([size(dw,1), 1]);
   else
      dw = [];
      ow = [];
   end
   actB = CS.activeCellFaces;
   actF = CS.activeFaces;
   nf = numel(orient) + numel(ow);

   Do = spdiags([orient; ow], 0, nf, nf) * blkdiag(CS.D(actB, actF), dw);
end

%--------------------------------------------------------------------------

function [totmob, omega] = dynamic_quantities(state, fluid)
   [mu, rho] = fluid.properties(state);
   s         = fluid.saturation(state);
   kr        = fluid.relperm(s, state);

   mob    = bsxfun(@rdivide, kr, mu);
   totmob = sum(mob, 2);
   omega  = sum(bsxfun(@times, mob, rho), 2) ./ totmob;
end
