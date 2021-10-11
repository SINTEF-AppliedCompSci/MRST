function state = solveIncompFlowLocal(state, g, s, fluid, varargin)
% Local version of solveIncompFlow for use with adjoint module. Local
% version has the additional option to supply right-hand-side to system
% directly. In addition, there is som extra slack in the assertion statement
% that checks sum(rates)=0  when only rate-controlled wells are present.
% This to prevent crash if input-constraint handling is not set to machine
% precision.
%
% SYNOPSIS:
%   state = solveIncompFlow(state, G, S, fluid)
%   state = solveIncompFlow(state, G, S, fluid, 'pn1', pv1, ...)
%
% DESCRIPTION:
%   This function assembles and solves a (block) system of linear equations
%   defining interface fluxes and cell and interface pressures at the next
%   time step in a sequential splitting scheme for the reservoir simulation
%   problem defined by Darcy's law and a given set of external influences
%   (wells, sources, and boundary conditions).
%
% REQUIRED PARAMETERS:
%   state  - Reservoir and well solution structure either properly
%            initialized from function 'initState', or the results from a
%            previous call to function 'solveIncompFlow' and, possibly, a
%            transport solver such as function 'explicitTransport'.
%
%   G, S   - Grid and (mimetic) linear system data structures as defined by
%            function 'computeMimeticIP'.
%
%   fluid  - Fluid data structure as described by 'fluid_structure'.
%
% OPTIONAL PARAMETERS:
%   wells  - Well structure as defined by function 'addWell'.  May be empty
%            (i.e., W = []) which is interpreted as a model without any
%            wells.
%
%   bc     - Boundary condition structure as defined by function 'addBC'.
%            This structure accounts for all external boundary conditions
%            to the reservoir flow.  May be empty (i.e., bc = []) which is
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
%   Solver - Which solver mode function 'solveIncompFlow' should employ in
%            assembling and solving the block system of linear equations.
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
%                              - pressure -- Well bottom-hole pressure.
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
%   `initState`, `solveIncompFlowMS`.

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


   opt = struct('bc', [], 'src', [], 'wells', [], 'rhs', [], ...
                'bcp',[],                                   ...
                'Solver',       'hybrid',                    ...
                'LinSolve',     @mldivide,                   ...
                'BlkDiag',      @blkdiag,                    ...
                'MatrixOutput', false,                       ...
                'pc_form','nonwetting');
   opt = merge_options(opt, varargin{:});

   % Check opt.Solver
   assert(sum(strcmpi(opt.Solver,  {'mixed', 'hybrid','tpfa'}))==1, ...
      ['Solver option must be one of ''tpfa'', ''mixed'', or ',...
       '''hybrid''.'])

   do_solve = check_input(g, opt);

   if do_solve,
      cellNo  = rldecode(1:g.cells.num, diff(g.cells.facePos), 2) .';
      s.C     = sparse(1:numel(cellNo), cellNo, 1);
      s.D     = sparse(1:numel(cellNo), double(g.cells.faces(:,1)), 1, ...
                       numel(cellNo), g.faces.num);
      s.sizeB = repmat(size(g.cells.faces, 1), [1,2]);
      s.sizeC = size(s.C);
      s.sizeD = size(s.D);

      [A, b, dF, dC] = build_system(state, g, s, opt.wells, ...
                                    opt.bc, opt.src, fluid, opt);

      solver   = pick_solver(g, s, dF, dC, opt);
      x        = solver(A, b);

      state = pack_solution(state, g, s, x{1:3}, opt);

      if opt.MatrixOutput, state.A = x{4}; end
   end
end

%--------------------------------------------------------------------------
% Helpers follow.
%--------------------------------------------------------------------------

function do_solve = check_input(g, opt)
   pressure_bc = ~isempty(opt.bc) && ...
                 any(strcmpi('pressure', opt.bc.type));
   well_bhp    = ~isempty(opt.wells) && ...
                 any(strcmpi('bhp', { opt.wells.type }));
   if ~(pressure_bc || well_bhp),
      sum_rates = 0;
      if ~isempty(opt.wells), sum_rates = sum(vertcat(opt.wells.val)); end
      if ~isempty(opt.src), sum_rates = sum_rates + sum(opt.src.rate);  end
      if ~isempty(opt.bc), sum_rates = sum_rates + sum(opt.bc.value);  end
      if abs(sum_rates) > 1e-5/day,
         error(id('MassBalance:NotFulfilled'), ...
              ['Well rates and flux bc must sum up to 0 \n', ...
               'when there are no bhp constrained wells or pressure bc.\n']);
      end
   end

   g_vec = gravity();
   % Check if there are gravity components in the grid axes or if the
   % pressure functions have been overridden (which means we cannot assume
   % anything about the influence of gravity and we play it safe).
   grav  = norm(g_vec(1 : g.griddim)) > 0 || isfield(g, 'grav_pressure');

   % We assemble and solve a system if there are any external forces or
   % sources (e.g., gravity, boundary conditions, explicit source terms,
   % injection/production wells or, in the case of the adjoint method, if
   % the caller supplied the right hand side directly).
   do_solve = grav || ~all([isempty(opt.bc),    ...
                            isempty(opt.src),   ...
                            isempty(opt.wells), ...
                            isempty(opt.bcp),   ...
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

function s = id(s)
   s = ['solveIncompFlow:', s];
end

%--------------------------------------------------------------------------

function [A, b, dF, dC] = build_system(state, g, s, w, bc, src, fluid, opt)
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
   %    A{1} = B = [ Br, 0; 0, Bw ]  % BLKDIAG
   %    A{2} = C = [ Cr   ;    Cw ]  % VERTCAT
   %    A{3} = D = [ Dr, 0; 0, Dw ]  % BLKDIAG
   %
   %    b{1} = f = [ fr   ;    fw ]  % VERTCAT
   %    b{2} = g = [ gr   ;       ]  % No direct contr. to cells from wells
   %    b{3} = h = [ hr   ;    hw ]  % VERTCAT
   %
   [Ar, br, dFr, dCr] = syscomp_res  (g, s, totmob, omega, bc, src, opt);
   %% made to add capillary pressure
   if(isfield(fluid,'pc'))
      pc=fluid.pc(state);
      if(~all(pc==0))
         mu  = fluid.properties(state);
         s   = fluid.saturation(state);
         kr  = fluid.relperm(s, state);
         mob = bsxfun(@rdivide, kr, mu);
         cc  = capPressureRHS(g, mob, pc, opt.pc_form);
         br{1} = br{1} + cc;
      end
   end

   [Aw, bw, dFw, dCw] = syscomp_wells(g, w, totmob, omega,          opt);

   A    = cell([1, 3]);               b    = cell([1, 3]);
   A{1} = opt.BlkDiag(Ar{1}, Aw{1});  b{1} = vertcat(br{1}, bw{1});
   A{2} = vertcat    (Ar{2}, Aw{2});  b{2} =         br{2}        ;
   A{3} = opt.BlkDiag(Ar{3}, Aw{3});  b{3} = vertcat(br{3}, bw{2});

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

function solver = pick_solver(g, s, dF, dC, opt)
   regul = ~any(dF);  % Set zero level if no prescribed pressure values.

   % Specify prescribed pressures on 'Dirichlet' contacts.
   % Initialize facePressure as NaN since the variable is only computed for
   % the boundary in the mixed solver.
   lam     = nan(size(dF));
   lam(dF) = dC;

   switch lower(opt.Solver),
      case 'hybrid',
         if ~any(strcmp(s.type, {'hybrid', 'comp_hybrid'})),
            error(id('SolverMode:Inconsistent'), ...
                 ['Solver mode ''hybrid'' is incompatible with ', ...
                  'linear system type ''%s''.'], s.type);
         end

         solver = @(A,b) solve_hybrid(A, b, lam, dF, regul, opt);
      case 'mixed',
         if ~any(strcmp(s.type, {'mixed', 'tpfa', 'comp_hybrid'})),
            error(id('SolverMode:Inconsistent'), ...
                 ['Solver mode ''mixed'' is incompatible with ', ...
                  'linear system type ''%s''.'], s.type);
         end

         solver = @(A,b) solve_mixed(A, b, g, s, lam, dF, ...
                                     regul, opt, @mixedSymm);
      case 'tpfa',
         if ~any(strcmp(s.type, {'mixed', 'tpfa', 'comp_hybrid'})),
            error(id('SolverMode:Inconsistent'), ...
                 ['Solver mode ''tpfa'' is incompatible with ', ...
                  'linear system type ''%s''.'], s.type);
         end

         solver = @(A,b) solve_mixed(A, b, g, s, lam, dF, ...
                                     regul, opt, @tpfSymm);
      otherwise,
         error(id('SolverMode:NotSupported'), ...
               'Solver mode ''%s'' is not supported.', opt.solver);
   end
end

%--------------------------------------------------------------------------

function x = solve_hybrid(A, b, lam, dF, regul, opt)
   A{3} = A{3}(:,~dF);
   b{3} = b{3}(  ~dF);
   [x{1:4}] = schurComplementSymm(A{:}, b{:},          ...
                                  'Regularize', regul, ...
                                  'LinSolve', opt.LinSolve);
   lam(~dF) = x{3};
   x{3}     = lam;
end

%--------------------------------------------------------------------------

function x = solve_mixed(A, b, g, s, lam, dF, regul, opt, solver)
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

%--------------------------------------------------------------------------

function state = pack_solution(state, g, s, flux, pres, lam, opt)
   if ~isfield(state, 'facePressure'),
      state.facePressure = zeros([g.faces.num, 1]);
   end
   state.pressure(:)     = pres;
   state.facePressure(:) = lam (1 : s.sizeD(2));
   state.flux(:)         = cellFlux2faceFlux(g, flux(1 : s.sizeB(1)));

   if ~isempty(opt.wells),
      nw  = numel(opt.wells);
      i_f = s.sizeB(1);
      i_p = s.sizeD(2);

      for k = 1 : nw,
         nperf = numel(opt.wells(k).cells);

         state.wellSol(k).flux     = -flux(i_f + 1 : i_f + nperf);
         state.wellSol(k).pressure =  lam (i_p + 1 : i_p +   1  );

         i_f = i_f + nperf;
         i_p = i_p +   1  ;
      end
   end
end

%--------------------------------------------------------------------------

function [A, b, dF, dC] = syscomp_res(g, s, mob, omega, bc, src, opt)
   A = cell([1, 3]);

   mob = spdiags(s.C * mob, 0, s.sizeB(1), s.sizeB(2));

   if strcmpi(opt.Solver, 'hybrid'),
      if ~isfield(s, 'BI'),
         error(id('SolverMode:Inconsistent'), ...
              ['Solver mode ''hybrid'' is incompatible with ', ...
               'linear system type ''%s''.'], s.type);
      end
      A{1} = mob * s.BI;
   else
      if ~isfield(s, 'B'),
         error(id('SolverMode:Inconsistent'), ...
              ['Solver mode ''%s'' is incompatible with ', ...
               'linear system type ''%s''.'], opt.Solver, s.type);
      end
      A{1} = mob \ s.B;
   end

   A{2} = s.C;
   A{3} = s.D;

   [b{1:3}, grav, dF, dC] = computePressureRHS(g, omega, bc, src);


   if ~isempty(opt.bcp)
      dp_face=zeros(g.faces.num,1);
      dp_face(opt.bcp.face)=opt.bcp.value;

      % use down stream pressure as variable, find the other
      % by global operations, assume the two corresponding
      % faces has different cartesian tag with same direction
      i=repmat([true,false],1,3);
      ind=find(i(g.cells.faces(:,2)));
      %dp_hfaces=zeros(size(g.cells.faces,1),1);
      dp_hfaces=dp_face(g.cells.faces(:,1))/2;
      dp_hfaces(ind)=-dp_hfaces(ind);
      b{1}=b{1}+dp_hfaces;
   end


   b{1} = b{1} + grav;  % Add (mimetic) gravity contributions.
end

%--------------------------------------------------------------------------

function [A, b, dF, dC] = syscomp_wells(G, W, mob, omega, opt)
   % Assume empty well structure by default...
   A  = cell([1, 3]);
   b  = cell([1, 2]);
   dF = logical([]);
   dC = [];

   if ~isempty(W),
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
      dp   = norm(gravity()) * vertcat(W.dZ) .* omega(wc);
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

   if ~isempty(opt.wells),
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
   if ~isempty(opt.wells),
      n  = sum(cellfun(@numel, { opt.wells.cells }));
      dw = speye(n);
      ow = ones([n, 1]);
   else
      dw = [];
      ow = sparse(0,0);
   end
   nf = numel(orient) + numel(ow);
   Do = spdiags([orient; ow], 0, nf, nf) * opt.BlkDiag(s.D, dw);
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
