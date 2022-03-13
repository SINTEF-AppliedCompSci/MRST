function state = solveIncompFlowVE(state, g, s, rock, fluid, varargin)
%Solve incompressible flow problem (fluxes/pressures) for VE equation.
%
% SYNOPSIS:
%   state = solveIncompFlowVE(state, G, S, rock, fluid)
%   state = solveIncompFlowVE(state, G, S, rock, fluid, 'pn1', pv1, ...)
%
% DESCRIPTION:
%   This function assembles and solves a (block) system of linear equations
%   defining interface fluxes and cell and interface pressures at the next
%   time step in a sequential splitting scheme for the reservoir simulation
%   problem defined by Darcy's law and the Vertical Equilibrium assumtion
%   and a given set of external influences (wells, sources, and boundary
%   conditions).
%
% NOTE:
%
%
% REQUIRED PARAMETERS:
%   state  - Reservoir and well solution structure either properly
%            initialized from functions 'initResSol' and 'initWellSol'
%            respectively, or the results from a previous call to function
%            'solveIncompFlowVE' and, possibly, a transport solver such as
%            function 'explicitTransportVE'.
%
%   G      - Grid as defined by function 'topSurfaceGrid'.
%
%   S      - (mimetic) linear system data structures as defined by
%            function 'computeMimeticIPVE'. NB: If system should be solved
%            with pseudo mobilities computed with vertical integration (and
%            not averaged permeabilites, then S must have been computed
%            using option 'intVert', true.
%
%   rock   - Rock data structure with valid field 'perm' from 3D model.
%
%   fluid  - Fluid object as defined by function 'initVEFluid'.
%
% OPTIONAL PARAMETERS:
%   wells  - Well structure as defined by function 'addWell'.  May be empty
%            (i.e., W = []) which is interpreted as a model without any
%            wells.
%
%   bc     - Boundary condition structure as defined by function 'addBC'.
%            This structure accounts for all external boundary conditions to
%            the reservoir flow.  May be empty (i.e., bc = []) which is
%            interpreted as all external no-flow (homogeneous Neumann)
%            conditions.
%
%   src    - Explicit source contributions as defined by function
%            'addSource'.  May be empty (i.e., src = []) which is
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
%   Solver - Which solver mode function 'solveIncompFlowVE' should employ
%            in assembling and solving the block system of linear equations.
%            String.  Default value: Solver = 'hybrid'.
%
%            Supported values are:
%              - 'hybrid' --
%                   Assemble and solve hybrid system for interface
%                   pressures.  System is eventually solved by Schur
%                   complement reduction and back substitution.
%
%                   The system 'S' must in this case be assembled by
%                   passing option pair ('Type','hybrid') or option pair
%                   ('Type','comp_hybrid') to function 'computeMimeticIP'.
%
%              - 'mixed' --
%                   Assemble and solve a hybrid system for interface
%                   pressures, cell pressures and interface fluxes. System
%                   is eventually reduced to a mixed system as per function
%                   'mixedSymm'.
%
%                   The system 'S' must in this case be assembled by
%                   passing option pair ('Type','mixed') or option pair
%                   ('Type','comp_hybrid') to function 'computeMimeticIP'.
%
%              - 'tpfa' --
%                   Assemble and solve a cell-centred system for cell
%                   pressures.  Interface fluxes recovered through back
%                   substitution.
%
%                   The system 'S' must in this case be assembled by
%                   passing option pair ('Type','mixed') or option pair
%                   ('Type','comp_hybrid') to function 'computeMimeticIP'.
%
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
%            caller of function 'solveIncompFlow'.
%            Logical.  Default value: MatrixOutput = FALSE.
%
% RETURNS:
%   state - Update reservoir and well solution structure with new values
%           for the fields:
%              - pressure -- Pressure values for all cells in the
%                            discretised reservoir model, 'G'.
%              - facePressure --
%                            Pressure values for all interfaces in the
%                            discretised reservoir model, 'G'.
%              - flux     -- Flux across global interfaces corresponding to
%                            the rows of 'G.faces.neighbors'.
%              - A        -- System matrix.  Only returned if specifically
%                            requested by setting option 'MatrixOutput'.
%
%              - wellSol  -- Well solution structure array, one element for
%                            each well in the model, with new values for
%                            the fields:
%                              - flux     -- Perforation fluxes through all
%                                            perforations for corresponding
%                                            well.  The fluxes are
%                                            interpreted as injection
%                                            fluxes, meaning positive
%                                            values correspond to injection
%                                            into reservoir while negative
%                                            values mean
%                                            production/extraction out of
%                                            reservoir.
%                              - bhp      -- Well bottom-hole pressure.
%
% NOTE:
%   If there are no external influences, i.e., if all of the structures
%   'W', 'bc', and 'src' are empty and there are no effects of gravity, and
%   no system right hand side has been supplied externally, then the input
%   values 'xr' and 'xw' are returned unchanged and a warning is printed in
%   the command window.  This warning is printed with message ID
%
%           'solveIncompFlow:DrivingForce:Missing'
%
% SEE ALSO:
%   `computeMimeticIP`, `addBC`, `addSource`, `addWell`, `initSimpleFluid`
%   `initResSol`, `initWellSol`, `solveIncompFlowMS`.

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

% $Date: 2012-09-06 09:34:34 +0200 (Thu, 06 Sep 2012) $
% $Revision: 9555 $

   opt = struct('bc', [], 'src', [], 'wells', [], 'rhs', [], ...
                'Solver',       'hybrid',                    ...
                'LinSolve',     @mldivide,                   ...
                'MatrixOutput', false);
   opt = merge_options(opt, varargin{:});

   do_solve = check_input(opt);

   if do_solve
      cellNo  = rldecode(1:g.cells.num, diff(g.cells.facePos), 2) .';
      s.C     = sparse(1:numel(cellNo), cellNo, 1);
      s.D     = sparse(1:numel(cellNo), double(g.cells.faces(:,1)), 1, ...
                       numel(cellNo), g.faces.num);
      s.sizeB = repmat(size(g.cells.faces, 1), [1,2]);
      s.sizeC = size(s.C);
      s.sizeD = size(s.D);

      [A, b, dF, dC] = build_system(state, g, s, opt.wells, ...
                                    opt.bc, opt.src, fluid, rock, opt);

      solver   = pick_solver(g, s, dF, dC, opt);
      x        = solver(A, b);

      state = pack_solution(state, g, s, x{1:3}, opt);
      if opt.MatrixOutput, state.A = x{4}; end
   else
      % No driving forces, set flux to zero and remove pressure
      state.flux(:) = 0;
      state.pressure = [];
      if isfield(state, 'facePressure')
         state.facePressure = [];
      end
   end
end

%--------------------------------------------------------------------------
% Helpers follow.
%--------------------------------------------------------------------------

function do_solve = check_input(opt)
   pressure_bc = ~isempty(opt.bc) && ...
                 any(strcmpi('pressure', opt.bc.type));
   well_bhp    = ~isempty(opt.wells) && ...
                 any(strcmpi('bhp', { opt.wells.type }));
   if ~(pressure_bc || well_bhp)
      sum_rates = 0;
      if ~isempty(opt.wells), sum_rates = sum(vertcat(opt.wells.val)); end
      if ~isempty(opt.bc), sum_rates = sum_rates + sum(opt.bc.value);  end
      if abs(sum_rates) > eps
         error(id('MassBalance:NotFulfilled'), ...
              ['Well rates and flux bc must sum up to 0 \n', ...
               'when there are no bhp constrained wells or pressure bc.\n']);
      end
   end

   grav  = norm(gravity) > 0;

   % We assemble and solve a system if there are any external forces or
   % sources (e.g., gravity, boundary conditions, explicit source terms,
   % injection/production wells or, in the case of the adjoint method, if
   % the caller supplied the right hand side directly).
   do_solve = grav || ~all([isempty(opt.bc),    ...
                            isempty(opt.src),   ...
                            isempty(opt.wells), ...
                            isempty(opt.rhs)]);
   if ~do_solve
      warning(id('DrivingForce:Missing'),                      ...
             ['\nNo external driving forces present in model--', ...
              'state remains unchanged: state.flux = 0, state.pressure = [].\n']);
   end

   % Assemble (and solve) system even in absence of external driving forces
   % if the caller requested 'MatrixOutput'.
   do_solve = do_solve || opt.MatrixOutput;
end

%--------------------------------------------------------------------------

function s = id(s)
   s = ['solveIncompFlow:', s];
end

%--------------------------------------------------------------------------

function [A, b, dF, dC] = build_system(state, g, s, w, bc, src, fluid, rock, opt)

   if (isfield(s, 'intVert') && s.intVert )
      % compute mobilites by vertically integrating 3D model
      % NB: should move computation of kr_H outside
      kr_H = integrateVertically(rock.perm(:,1), inf, g);
      kr = integrateVertically(rock.perm(:,1), state.h, g);
      mob = [kr./fluid.mu(1), (kr_H-kr)./fluid.mu(2)];
   else
      mob=  fluid.mob(state);
   end

   totmob   = sum(mob,2);
   omega = sum(bsxfun(@times,mob,fluid.rho),2)./totmob;

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
   %    A{1} = B = [ Br, 0; 0, Bw ]  % BLKDIAG
   %    A{2} = C = [ Cr   ;    Cw ]  % VERTCAT
   %    A{3} = D = [ Dr, 0; 0, Dw ]  % BLKDIAG
   %
   %    b{1} = f = [ fr   ;    fw ]  % VERTCAT
   %    b{2} = g = [ gr   ;       ]  % No direct contr. to cells from wells
   %    b{3} = h = [ hr   ;    hw ]  % VERTCAT
   %
   [Ar, br, dFr, dCr] = syscomp_res  (g, s, totmob, omega, fluid.pc, state, bc, src, opt);
   [Aw, bw, dFw, dCw] = syscomp_wells(g, w, totmob, omega, fluid.pc, state, opt);

   A    = cell([1, 3]);           b    = cell([1, 3]);
   A{1} = blkdiag(Ar{1}, Aw{1});  b{1} = vertcat(br{1}, bw{1});
   A{2} = vertcat(Ar{2}, Aw{2});  b{2} =         br{2}        ;
   A{3} = blkdiag(Ar{3}, Aw{3});  b{3} = vertcat(br{3}, bw{2});

   % If rhs is supplied by user, replace b constructed above by given rhs
   if ~isempty(opt.rhs)
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

function solver = pick_solver(g, s, dF, dC, opt)
   regul = ~any(dF);  % Set zero level if no prescribed pressure values.

   switch lower(opt.Solver)
      case 'hybrid'
         if ~any(strcmp(s.type, {'hybrid', 'comp_hybrid'}))
            error(id('SolverMode:Inconsistent'), ...
                 ['Solver mode ''hybrid'' is incompatible with ', ...
                  'linear system type ''%s''.'], s.type);
         end
         
         solver = @solve_hybrid;
%       case 'mixed',
%          if ~any(strcmp(s.type, {'mixed', 'tpfa', 'comp_hybrid'})),
%             error(id('SolverMode:Inconsistent'), ...
%                  ['Solver mode ''mixed'' is incompatible with ', ...
%                   'linear system type ''%s''.'], s.type);
%          end
%
%          solver = @(A,b) solve_mixed(A, b, @mixedSymm);
%       case 'tpfa',
%          if ~any(strcmp(s.type, {'mixed', 'tpfa', 'comp_hybrid'})),
%             error(id('SolverMode:Inconsistent'), ...
%                  ['Solver mode ''tpfa'' is incompatible with ', ...
%                   'linear system type ''%s''.'], s.type);
%          end
%
%          solver = @(A,b) solve_mixed(A, b, @tpfSymm);
       otherwise
         error(id('SolverMode:NotSupported'), ...
               'Solver mode ''%s'' is not supported.', opt.solver);
   end

   % Specify prescribed pressures on 'Dirichlet' contacts.
   % Initialize facePressure as NaN since the variable is only computed for
   % the boundary in the mixed solver.
   lam     = nan(size(dF)); %zeros(size(dF));
   lam(dF) = dC;

   function x = solve_hybrid(A, b)
      A{3} = A{3}(:,~dF);
      b{3} = b{3}(  ~dF);
      [x{1:4}] = schurComplementSymm(A{:}, b{:},          ...
                                     'Regularize', regul, ...
                                     'LinSolve', opt.LinSolve);
      lam(~dF) = x{3};
      x{3}     = lam;
   end

   function x = solve_mixed(A, b, solver)
      nF = neumann_faces(g, dF, opt);
      cellNo = rldecode(1 : g.cells.num, diff(g.cells.facePos), 2) .';
      sgn = 2*double(g.faces.neighbors(g.cells.faces(:,1), 1) == cellNo) - 1;

      Do = oriented_mapping(s, sgn, opt);

      A{3} = A{3}(:,nF);
      b{3} = b{3}(  nF);
      [x{1:4}] = solver(A{:}, b{:}, Do, 'Regularize', regul, ...
                        'LinSolve', opt.LinSolve);

      x{1}    = [faceFlux2cellFlux(g, x{1}(1 : g.faces.num)); ...
                 x{1}(g.faces.num + 1 : end)];

      lam(nF) = x{3};
      x{3}    = lam;
   end
end

%--------------------------------------------------------------------------

function state = pack_solution(state, g, s, flux, pres, lam, opt)
   if ~isfield(state, 'facePressure')
      state.facePressure = zeros([g.faces.num, 1]);
   end
   state.pressure(:)     = pres;
   state.facePressure(:) = lam (1 : s.sizeD(2));
   state.flux(:)         = cellFlux2faceFlux(g, flux(1 : s.sizeB(1)));

   if ~isempty(opt.wells)
      nw  = numel(opt.wells);
      i_f = s.sizeB(1);
      i_p = s.sizeD(2);

      for k = 1 : nw
         nperf = numel(opt.wells(k).cells);

         state.wellSol(k).flux     = -flux(i_f + 1 : i_f + nperf);
         state.wellSol(k).bhp =  lam (i_p + 1 : i_p +   1  );

         i_f = i_f + nperf;
         i_p = i_p +   1  ;
      end
   end
end

%--------------------------------------------------------------------------

function [A, b, dF, dC] = syscomp_res(g, s, mob, omega, pc, state, bc, src, opt)
   A = cell([1, 3]);

   mob = spdiags(s.C * mob, 0, s.sizeB(1), s.sizeB(2));

   if strcmpi(opt.Solver, 'hybrid')
      if ~isfield(s, 'BI')
         error(id('SolverMode:Inconsistent'), ...
              ['Solver mode ''hybrid'' is incompatible with ', ...
               'linear system type ''%s''.'], s.type);
      end
      A{1} = mob * s.BI;
   else
      if ~isfield(s, 'B')
         error(id('SolverMode:Inconsistent'), ...
              ['Solver mode ''%s'' is incompatible with ', ...
               'linear system type ''%s''.'], opt.Solver, s.type);
      end
      A{1} = mob \ s.B;
   end

   A{2} = s.C;
   A{3} = s.D;

   [b{1:3}, grav, dF, dC] = computePressureRHSVE(g, omega, pc, bc, src, state);
   b{1} = b{1} + grav;  % Add (mimetic) gravity contributions.
end

%--------------------------------------------------------------------------

function [A, b, dF, dC] = syscomp_wells(G, W, mob, omega,pc, state, opt)
   % Assume empty well structure by default...
   A  = cell([1, 3]);
   b  = cell([1, 2]);
   dF = logical([]);
   dC = [];

   if ~isempty(W)
      % but fill in values when there nevertheless are some...
      nperf = cellfun(@numel, { W.cells });
      wc    = vertcat(W.cells);   n = numel(wc);   i = 1 : n;
      nW    = numel(W);

      % Diagonal transmissibility matrix (inner product).
      v = mob(wc) .* vertcat(W.WI);
      if ~strcmpi(opt.Solver, 'hybrid'), v = 1 ./ v; end
      A{1} = sparse(i, i, v, n, n); % == spdiags(v,0,n,n), but more direct.

      % Connection matrices {2} -> C and {3} -> D for all wells.
      A{2} = sparse(i, wc                      , 1, n, G.cells.num);
      A{3} = sparse(i, rldecode(1:nW, nperf, 2), 1, n, nW         );

      % Which (and what) are the prescribed well bottom-hole pressures?
      dF   = strcmpi('bhp', { W.type } .');
      dC   = reshape([ W(dF).val ], [], 1);

      % Form linsys rhs contributions, {1} -> pressure, {2} -> rate.
      b{1} = -vertcat(W.val);  b{1}(~dF) = 0;  % Remove rates
      b{2} = -vertcat(W.val);  b{2}( dF) = 0;  % Remove pressures

      % Expand well pressure rhs to each perforation, adjust for gravity
      % effects if applicable (i.e., when NORM(gravity())>0 and W.dZ~=0).
      pc_cells=pc(state.h);
      dp   = norm(gravity()) * ((vertcat(W.dZ) +pc_cells(wc)).* omega(wc));
      b{1} = rldecode(b{1}, nperf) - dp;
   end
end

%--------------------------------------------------------------------------

function nF = neumann_faces(g, dF, opt)
% Determine the 'Neumann faces' (nF) (i.e., the faces for which
% face/contact pressures will be determined in the 'mixed' and 'TPFA'
% cases) of the model (g) according to the following rules:
%
%   - No internal faces are counted amongst the 'Neumann faces'.
%   - An external face is a 'Neumann face' unless there is a prescribed
%     pressure value (i.e., a Dirichlet condition) associated to the face.
%   - Wells controlled by rate-constraints are treated as Neumann contacts.
%
   nF                                 = false([g.faces.num, 1]); % No int.
   nF(any(g.faces.neighbors == 0, 2)) = true;  % All external faces ...
   nF(dF(1 : g.faces.num))            = false; % ... except Dirichlet cond.

   if ~isempty(opt.wells)
      % Additionally include all rate-constrained wells.
      nF = [nF; strcmpi('rate', { opt.wells.type } .')];
   end
end

%--------------------------------------------------------------------------

function Do = oriented_mapping(s, orient, opt)
% 'Do' maps face fluxes to half-face fluxes.  This matrix is used to form
% the reduced, mixed system of linear equtions in linear system solver
% functions 'mixedSymm' and 'tpfSymm'.
%
   if ~isempty(opt.wells)
      n  = sum(cellfun(@numel, { opt.wells.cells }));
      dw = speye(n);
      ow = ones([n, 1]);
   else
      dw = [];
      ow = [];
   end
   nf = numel(orient) + numel(ow);
   Do = spdiags([orient; ow], 0, nf, nf) * blkdiag(s.D, dw);
end
