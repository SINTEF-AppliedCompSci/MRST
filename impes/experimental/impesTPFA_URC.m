function [state,dt,report,sreport] = impesTPFA_URC(state0, G, T, fluid, dt, pv, varargin)
%Single IMPES step based on approximate Newton pressure solve
%
% SYNOPSIS:
%   state = impesTPFA(state, G, T, fluid, dt, pv)
%   state = impesTPFA(state, G, T, fluid, dt, pv, 'pn1', pv1, ...)
%
% DESCRIPTION:
%   This function assembles and solves a non-linear system of equation that
%   defines phase pressures under the assumption of mass balance at some
%   future time.  The formulation uses the Newton (i.e., a residual and its
%   associated (approximate) Jacobian matrix) formalism in the system
%   definition.
%
%   Once the non-linear system has been solved to convergence, fluids are
%   displaced, using an explicit Euler method, according to the new mass
%   balance relations.
%
%   This function uses a two-point flux approximation (TPFA) method with
%   minimal memory consumption within the constraints of operating on a
%   fully unstructured polyhedral grid structure.  Moreover, we attempt to
%   minimise the number of fluid matrix evaluations.
%
% REQUIRED PARAMETERS:
%   state - Reservoir and well solution structure either properly
%           initialized from function 'initResSolComp', or the results from
%           a previous call to function 'compTPFA' and, possibly, a
%           transport solver such as function 'implicitTransport'.
%
%   G, T  - Grid and half-transmissibilities as computed by the function
%           'computeTrans'.
%
%   fluid - Compressible fluid object.
%
%   dt    - Time step size (IMPES/pressure step).
%
%   pv    - Pore volume.  One positive scalar for each grid cell.  Use
%           function 'poreVolume' to compute this value.
%
% OPTIONAL PARAMETERS (supplied in 'key'/value pairs ('pn'/pv ...)):
%   wells   - Well structure as defined by functions 'addWell'.  May be
%             empty (i.e., W = struct([])) which is interpreted as a model
%             without any wells.
%
%   bc      - Boundary condition structure as defined by function 'addBC'.
%             This structure accounts for all external boundary conditions
%             to the reservoir flow.  May be empty (i.e., bc = struct([]))
%             which is interpreted as all external no-flow (homogeneous
%             Neumann) conditions.
%
%   src     - Explicit source contributions as defined by function
%             'addSource'.  May be empty (i.e., src = struct([])) which is
%             interpreted as a reservoir model without explicit sources.
%
%   Verbose - Enable output.  Default value dependent upon global verbose
%             settings of function 'mrstVerbose'.
%
%   ATol, RTol -
%             Non-linear iteration convergence controls.  The computational
%             process is terminated when the infinity norm of the
%             non-linear residual is less than
%
%                      MAX(ATol, RTol*NORM(Initial Residual, INF))
%
%             Positive scalars.  Default values: ATol=5.0e-7, RTol=1.0e-8.
%
%   LinSolve -
%             Handle to linear system solver software to which the fully
%             assembled pressure updated system of linear equations will be
%             passed on each iteration of the Newton-Raphson procedure.
%             Assumed to support the syntax
%
%                        x = LinSolve(A, b)
%
%             in order to solve a system Ax=b of linear equations.
%             Default value: LinSolve = @mldivide (backslash).
%
%   UpdateMass -
%             Whether or not to update masses at the end of the pressure
%             solve.  Logical.  Default value: UpdateMass=true.  Set
%             'UpdateMass' to false to use solver in a sequential splitting
%             framework.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Options only for the URC type transport
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   URC_p    - Use URC option for pressure. include only one newton
%              iteration, and elimnation of some terms in the jacobian
%              If URC_p=false and max_it==1 linearizing with all terms in
%              the jacobian is achived
%
%   URC_s    - Use URC timestep update by calculating saturation
%              implicite.
%              acceptable values {'none','linear','nonlinear'}
%              if 'linear' only one newton iteration is tried.
%              if 'nonlinear' a full nonlinear solve is tried. The options
%              for this solve should be given in URC_options
%
%  URC_options - options for a possible full nonlinear solve of saturation
%
%  ntol     - tolerance for the nonliner saturaiton update, for now also used
%              to accept timestep
%
%  volume_error - tolerance for acceptable volume error after surface
%                 volume update
%
%  pressure_change - acceptable maximum pressure change of cell pressures
%
%
% RETURNS:
%   state - Updated reservoir and well solution structure with new values
%           for the fields:
%              - pressure -- Pressure values for all cells in the
%                            discretised reservoir model, 'G'.
%              - facePressure --
%                            Pressure values for all interfaces in the
%                            discretised reservoir model, 'G'.
%              - flux     -- Flux across global interfaces corresponding to
%                            the rows of 'G.faces.neighbors'.
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
%              - z        -- Surface volumes per phase.  One column for
%                            each fluid phase, and one row for each grid
%                            cell.
%
%     dt       - the time step actually used or at least tried
%
%     report   - report of well data
%
%     sreport  - struct giving information of the solution, if
%                sreport.success=true the timestep succeded else ther
%                were som error. Intended to be used for time step control
%
%
%
% SEE ALSO:
%   compTPFA, poreVolume.

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

   opt = struct('bc', [], 'src', [], 'wells', [],             ...
                'ATol', 5.0e-7, 'RTol', 1.0e-8,               ...
                'max_it',10,...
                'UpdateMass', true, 'DynamicMobility', false, ...
                'EstimateTimeStep', true,                     ...
                'LinSolve', @mldivide, 'Verbose', mrstVerbose, ...
                'LineSearch', false,...
                'wellModel', @impesTPFADefaultWellModel, ...
                'report', [],...
                'init_state',state0,...
                'face_update',true,...
                'cfl_factor',1,...
                'URC_p', false,...%deside iterations and jacoby terms
                'URC_s', 'linear',...%deside if implicite transport
                'URC_opt',struct('maxnewt',10,'lstrials',3),...
                'ntol' , 1e-3,... %nonlinear error for transport
                'volume_error',1e-3,...%deside acceptable volume error
                'pressure_change',30*barsa,...% deside acceptable pressure change
                'rock',[]);%needed for implicite transport
   sreport=struct('success',true,...
      'dt',dt,...
      'p_E',0,...
      'upwind_changed',false,...
      'transport_failed',false,...
      'pressure_failed',false,...
      'max_volume_error',0,...
      'treport',[]);
   opt = merge_options(opt, varargin{:});
   state=opt.init_state;
   opt.wdp = opt.wellModel(state, opt.wells, fluid);

   trans         = compute_trans(G, T, opt);
   [cmob, cdmob] = impesComputeMobility(state, fluid, opt.bc, ...
                                        opt.wells, opt.wdp);

   OP = tpfaConnectionOperators(G, opt.wells, size(state.z, 2));

   [luAc, dAc, Af, mob, dmob, density] = ...
      eval_fluid_data(state, G, fluid, cmob, cdmob, opt, OP);

   save = struct('luAc', luAc, 'dAc', dAc, 'Af', Af, ...
                 'mob', mob, 'dmob', dmob);
   if opt.EstimateTimeStep,
      rho = []; pc = []; dpc = []; cfl_fac = 1;
      ix  = 1 : G.faces.num;

      est_dt = @(press, mob, dmob) ...
         opt.cfl_factor*estimate_dt_coats(G, trans(ix), pv, mob(ix, :), dmob(ix, :),  ...
                           press, rho, pc, dpc, cfl_fac);

      clear ix rho pc dpc cfl_fac
   end

   [F, fcontrib_dp, wcontrib_dp, gflux] = ...
      assemble_residual(luAc, Af, mob, state, state0,  G, trans, dt, pv, ...
                        density, opt, OP, save.luAc);
   E0   = norm(F, inf);   E = E0;
   done = E < opt.ATol;
   it   = 0;
   dispif(opt.Verbose && ~done, ...
         ['Solve pressure system by approximate Newton method.\n', ...
          'Terminate when NORM(F,INF) < MAX(%8.2e, %8.2e)\n'], ...
          opt.ATol, opt.RTol * E0);
   while ~done,
      if opt.EstimateTimeStep,
         dt_tmp = est_dt(state.pressure, mob, dmob);
         dt = min(dt, dt_tmp);
      end
      [F, fcontrib_dp, wcontrib_dp, gflux] = ...
      assemble_residual(luAc, Af, mob, state, state0,  G, trans, dt, pv, ...
                        density, opt, OP, save.luAc);
      E    = norm(F, inf);

      done = (E < opt.ATol) || (E < opt.RTol * E0) || ...
             (opt.URC_p && (it > 0) || it > opt.max_it);  % == opt.URC ...
      dispif(opt.Verbose, ...
             '|| F(%02d) ||_oo = %11.5e\n', it, E);
      if done
         continue;
      end

      J = approximate_jacobian(luAc, dAc, G, state, state0,     ...
                               fcontrib_dp, wcontrib_dp, ...
                               dt, pv, opt, OP, save.luAc, ...
                               mob, fluid);

      dpress = - opt.LinSolve(J, F);


      if opt.LineSearch,
         dpress = line_search(dpress, G, T, trans, dt, state, state0, fluid, ...
                              mob, density, cmob, cdmob, pv, F, opt, OP);
      end


      state = update_pressure(state, dpress, G, T, mob, density, OP, opt);

%{
      opt.wdp = opt.wellModel(state, opt.wells, fluid);
%}

      if opt.DynamicMobility,
         [cmob, cdmob] = impesComputeMobility(state, fluid, opt.bc, ...
                                              opt.wells, opt.wdp);
      end

      [luAc, dAc, Af, mob, dmob, density] = ...
         eval_fluid_data(state, G, fluid, cmob, cdmob, opt, OP);

      save.luAc = luAc;

      % luAc = save.luAc;
      % dAc  = save.dAc;
      % mob  = save.mob;
      if(~opt.face_update)
          Af   = save.Af;
      end

      it   = it + 1;
   end


   s      = { 'ATol', 'RTol' };
   plural = { ''    , 's'    };
   dispif(opt.Verbose, ...
         ['Pressure converged in %d iteration%s (%s criterion): ', ...
          '|| F(end) ||_oo = %12.5e\n'], it, plural{1 + (it ~= 1)}, ...
          s{1 + ~(E < opt.ATol)}, E);

   [state, dp, flux] = compute_flux(state, G, trans, density, mob, opt, OP);
   sreport.success=true;
   sreport.dt=dt;
   if(~strcmp(opt.URC_s,'none'))
      % check if resonable accuracy
      sreport.p_E=E;
      if(~(E <opt.ATol))
        sreport.success=false;
        sreport.pressure_failed=true;
        dispif(opt.Verbose,'Pressure step failed: to large error',E);
        if(opt.Verbose)
           disp(sreport);
        end
      end
      % check if pressure has changed to much
      opt.max_dp=max(abs(state.pressure-state.pressure>opt.pressure_change));
      if(opt.max_dp>opt.pressure_change)
         sreport.success=false;
         sreport.pressure_failed=true;
         dispif(opt.Verbose,'Pressure step failed : To large pressure change: ',max(abs(state.pressure-state.pressure)),'\n');
      end
      % calculate upwind directions from pressure solution
      % should assert that this is equivalent with phase pressure
      % difference
      ix_pressure=findUpWindFromFlux(G,state,fluid,opt.rock,opt);
      switch opt.URC_s
       case 'linear'
        [state_new, wreport, report]  = implicitTransportSatURC(state, G, dt, opt.rock, fluid, 'wells', opt.wells,...
                                                          'verbose' , mrstVerbose , ...  % Emit progress reports?
                                                          'nltol'   , opt.ntol, ...  % Non-linear residual tolerance
                                                          'lstrials', 1    , ...  % Max no of line search trials
                                                          'maxnewt' , 1    , ...  % Max no. of NR iterations
                                                          'force_iteration', true,...%always take a newton step
                                                          'resred'  , 1,...
                                                          'init_state', ...
                                                          state);
       case 'nonlinear'
        [state_new, wreport, report]  = implicitTransportSatURC(state, G, dt, opt.rock, fluid, 'wells', opt.wells,...
                                                          'verbose' , mrstVerbose , ...  % Emit progress reports?
                                                          'nltol'   , opt.ntol, ...  % Non-linear residual tolerance
                                                          'lstrials', opt.URC_opt.lstrials    , ...  % Max no of line search trials
                                                          'maxnewt' , opt.URC_opt.maxnewt    , ...  % Max no. of NR iterations
                                                          'force_iteration', true,...%always take a newton step
                                                          'resred'  , 1,...
                                                          'init_state', ...
                                                          state);
       otherwise
        error('No such URC method')
      end
        % calculate upwind directions after transport
       ix_transport=findUpWindFromFlux(G,state,fluid,opt.rock,opt);
       %     if upwind direction has changed,
       %        NB
       if(~all(ix_transport(:)==ix_pressure(:)))
          warning('Upwind directions change: URC step failed');
          sreport.success=false;
          sreport.upwind_change=true;
          % we should no do pressure once more with new saturations  to
          % evaluate mobility
          sreport.init_state=state_new;
          % here we potentialy use the new z values also for evaluating
          % propertice, but for now it is nesseary since code tend to use z
          % to calculate saturation
       end

        % we should check:
        %     solution is acceptable
        sreport.treport=report;
        if(report.success==false);
           sreport.success=false;
           sreport.transport_failed=true;
           dispif(opt.Verbose,'impliciteTransport failed\n')
        end

        % make a state with the new saturation but with old z values
        % mobility is evaluated in new mobilities but viscosity is taken
        % based on old z
        state.s=state_new.s;
        [cmob, cdmob] = impesComputeMobilityURC(state, fluid, opt.bc, ...
           opt.wells, opt.wdp);
        % use mobility from new saturation but old z for RS
        %[luAc, dAc, Af, mob, dmob, density] = ...
        [luAc, dAc, Af, mob] = ...
           eval_fluid_data(state, G, fluid, cmob, cdmob, opt, OP);
        [state_tmp, ok] = update_masses_flux(state0.z, flux, Af, mob, dt, pv, OP);
        if(ok==false)
           sreport.success=false;
           sreport.mass_update_failed=true;
           dispif(opt.Verbose,'Update masses failed\n')
        end
        % check if masses are ok
        state.z=state_tmp.z;
        [u, u, u, u] = fluid.pvt(state.pressure, state.z);
        % check if volume has small error
        sreport.max_volume_error=max(abs(sum(u,2)-1));
        if(sreport.max_volume_error>opt.volume_error)
           sreport.success=false;
           dispif(opt.Verbose,'To large volume errror\n')
        end
        state.s      = bsxfun(@rdivide, u, sum(u, 2));
        if(sreport.success==false)
           if(opt.Verbose)
              disp('impes URC step failed')
              disp(sreport)
           end
           return;
        end
        %assert(all(state.s==s0))
   elseif(opt.UpdateMass),
      state_tmp = update_masses(state0, dp, fcontrib_dp, gflux, dt, pv, OP);
      state.z=state_tmp.z;
      [u, u, u, u] = fluid.pvt(state.pressure, state.z);               %#ok
      state.s      = bsxfun(@rdivide, u, sum(u, 2));
      % update face mobility and Af for the definition of
      % well rates
      [cmob, cdmob] = impesComputeMobility(state, fluid, opt.bc, ...
      opt.wells, opt.wdp);
      %[luAc, dAc, Af, mob, dmob, density] = ...
      [luAc, dAc, Af, mob] = ...
         eval_fluid_data(state, G, fluid, cmob, cdmob, opt, OP);
   end
   % this is computed to make the report corrct
   assert(all(state.z(:)>=0))
   assert(all(state.pressure(:)>=0))
   if(~isempty(opt.report))
      time=dt+opt.report(end).TIME(end);
   else
      time=dt;
   end
   offset = G.faces.num*size(state.z, 2);
   report = [opt.report; ...
      makeReport(state, G, mob, Af(offset+1:end, offset+1:end), time, it)];
end

%--------------------------------------------------------------------------

function report = makeReport(state, G, mob, Aw, time, nit)
   if ~isfield(state, 'wellSol') || isempty(state.wellSol),
      report = [];
      return;
   end

   wno = rldecode(1 : numel(state.wellSol), ...
                  cellfun('prodofsize', { state.wellSol.flux }), 2) .';

   Accum = sparse(wno, 1 : numel(wno), 1);
   wmob  = mob((G.faces.num + 1) : end, :);
   f     = bsxfun(@rdivide, wmob, sum(wmob, 2));

   perfflux = -vertcat(state.wellSol.flux);
   resrate  = bsxfun(@times, f, perfflux);
   surfrate = mmultiply(Aw, resrate);

   arates = Accum * [surfrate, resrate, perfflux];

   TIME = time;
   WBHP = vertcat(state.wellSol.pressure);

   report = struct('TIME', TIME, 'NEWT', nit, 'WBHP', WBHP);

   if size(mob, 2) == 3,
      % Assume Aqua <-> 1, Liquid <-> 2, Vapour <-> 3
      %
      WWPR = arates(:, 1);
      WOPR = arates(:, 2);
      WGPR = arates(:, 3);

      WVPT = arates(:,  end );
      WWCT = arates(:, 3 + 1) ./ WVPT;
      WGCT = arates(:, 3 + 3) ./ WVPT;

      report.WVPT = WVPT;  report.WWPR = WWPR;  report.WOPR = WOPR;
      report.WGPR = WGPR;  report.WWCT = WWCT;  report.WGCT = WGCT;
   end
end

%--------------------------------------------------------------------------

function trans = compute_trans(G, T, opt)
   trans = 1 ./ accumarray(G.cells.faces(:,1), 1 ./ T, [G.faces.num, 1]);

   if ~isempty(opt.wells),
      trans = [trans; vertcat(opt.wells.WI)];
   end
end

%--------------------------------------------------------------------------

function [luAc, dAc, Af, mob, dmob, rho] = ...
      eval_fluid_data(state, G, fluid, mob, dmob, opt, OP)

   n = G.cells.num * size(mob, 2);

   [p, z, face_p] = impesAssembleStateVars(state, opt.bc, ...
                                           opt.wells, opt.wdp);
   [Ac, dAc, rho] = fluid_matrix(p, z, fluid);

   if size(Ac, 1) > n,
      % Possibly very costly.  Measure and investigate alternatives.

      i   = 1 : n;
      Ac  = Ac (i, i);
      dAc = dAc(i, i);
   end

   luAc = factorise(Ac);

   [mob, face_z, dmob] = tpfaUpwindStateVars(G, p, z, rho, mob, dmob, ...
                                             opt.bc, opt.wells, OP);

   Af = fluid_matrix(face_p, face_z, fluid);
end

%--------------------------------------------------------------------------

function [F, fcontrib_dp, wcontrib_dp, grav_term] = ...
      assemble_residual(luAc, Af, mob, state, state0, G, trans, ...
                        dt, pv, rho, opt, OP, luAc_tmp)

   dyntrans    = bsxfun(@times, trans, mob);
   fcontrib_dp = mmultiply(Af, dyntrans);
   grav_term   = grav_cap(Af, dyntrans, OP, rho(1:G.cells.num, :));
   press       = impesAssembleStateVars(state, opt.bc, opt.wells, opt.wdp);
   dp          = pressure_diff(G, press, opt);
   fcontrib    = bsxfun(@times, fcontrib_dp, dp) + grav_term;
   fcontrib(OP.connno(~OP.active),:)=0;

   F = residual_cell(luAc, state0, G, fcontrib, dt, pv, opt, OP, luAc_tmp);

   wcontrib_dp = [];

   if ~isempty(opt.wells),
      i = G.faces.num + 1 : size(dyntrans, 1);
      wdyntrans_r = dyntrans   (i, :);
      wdyntrans_s = fcontrib_dp(i, :);

      [Fw, wcontrib_dp] = residual_well(opt.wells, opt.wdp, state, ...
                                        wdyntrans_s, wdyntrans_r);
      %F  = [F; Fw];
      F  = [F; dt*Fw];
      wcontrib_dp=dt*wcontrib_dp;
   end
end

%--------------------------------------------------------------------------

function dp = line_search(dp0, G, T, trans, dt, state_old, state0, fluid, ...
                          mob, rho, cmob, cdmob, pv, F, opt, OP)
   norm_F0 = norm(F, inf);
   norm_F  = 10 * norm_F0;

   alpha   = 0;

   while ~(norm_F <   norm_F0),
      dp = pow2(dp0, alpha);
      dp = min(dp, 50*barsa);
      dp = max(dp, -50*barsa);
      state   = update_pressure(state_old, dp, G, T, mob, rho, OP, opt);

      [luAc, dAc, Af, mob, dmob, density] = ...
         eval_fluid_data(state, G, fluid, cmob, cdmob, opt, OP);  %#ok

      F = assemble_residual(luAc, Af, mob, state, state0, G, trans, dt, pv, ...
                            density, opt, OP, luAc);

      alpha = alpha - 1;
      norm_F = norm(F, inf);
   end

   dispif(mrstVerbose, 'alpha = %d\n', alpha + 1);
end

%--------------------------------------------------------------------------

function state = update_pressure(state, dp, G, T, mob, rho, OP, opt)
   cf  = G.cells.faces(:,1);
   p   = state.pressure + dp(1 : G.cells.num);
   dG  = bsxfun(@times, OP.gpot, rho(OP.cellno, :));

   num = bsxfun(@plus, p(OP.cellno), dG);
   Tm  = bsxfun(@times, T, mob(cf, :));
   fp  = accumarray(cf, sum(Tm .* num, 2)) ./ accumarray(cf, sum(Tm, 2));

   if ~isempty(opt.bc),
      di = strcmpi(opt.bc.type, 'pressure');
      fp(opt.bc.face(di)) = opt.bc.value(di);
   end

   i = find(~isfinite(fp));
   if ~isempty(i),
      % Assertion:
      %
      %   This situation occurs only in multi-phase flow and only if *all*
      %   (upwind) phase mobilites are identically zero on the faces 'i' in
      %   which case the above defintion produces fp=0/0 (==NaN).
      %   Moreover, these faces *must* be internal (not boundary).
      %
      %   A typical case in which this may happen is if the interfaces 'i'
      %   represent *sharp* separators between mobile phases of different
      %   mass densities with the heaviest at "the bottom".  Fairly
      %   academic, but nevertheless possible as a result of equilibrium
      %   based intialisation or as a stable state at T_\infty.
      %
      % Use a simple transmissibility weighted average of the cell
      % pressures in this case as there is effectively *no* connection
      % across the interface.  In other words, the exact interface pressure
      % value does not matter, we need only compute something vaguely
      % reasonable in order to subsequently compute A(fp).
      %
      % On the other hand, as there is no connection across the interface,
      % the Jacobian matrix will effectively decouple into sub-systems and
      % we may experience difficulty solving the system of linear eqns.
      %
      N = G.faces.neighbors(i, :);

      % Verify above assertion.
      %
      assert (size(mob, 2) > 1, ...
             ['Non-finite connection pressure in single-phase ', ...
              'flow is highly unexpected.']);

      assert (all(isnan(fp(i))), ...
              'Inifinite connection pressures highly unexpected.');

      assert (all(all(mob(i, :) == 0)), ...
             ['Non-finite connection pressure even if upwind ', ...
              'phase mobilities are non-zero.']);

      assert (all(all(N ~= 0, 2)), ...
              'Unexpected zero (upwind) total mobility on boundary.');

      % Compute connection pressure by means of transmissibility weighted
      % cell pressure values.
      %
      c  = reshape(N, [], 1);
      hf = OP.hf(repmat(i, [2, 1]), c);
      ix = repmat((1 : numel(i)).', [2, 1]);

      fp(i) = accumarray(ix, T(hf) .* p(c)) ./ accumarray(ix, T(hf));
   end

   if ~isempty(opt.wells),
      is_bhp = strcmpi('bhp', { opt.wells.type });

      pw          = zeros([numel(is_bhp), 1]);
      pw( is_bhp) = vertcat(opt.wells(is_bhp).val);

      if any(~is_bhp),
         pw(~is_bhp) = vertcat(state.wellSol(~is_bhp).pressure) + ...
                       dp(G.cells.num + 1 : end);
      end

      if ~all(pw > 0),
         error('Non-physical pressure increment in wells');
      end

      for w = 1 : numel(opt.wells),
         state.wellSol(w).pressure = pw(w);
      end
   end

   if ~all(p > 0),
      error('Non-physical pressure increment in cells');
   end

   state.pressure     = p;
   state.facePressure = fp;
end

%--------------------------------------------------------------------------

function [state, dp, flux] = compute_flux(state, G, trans, rho, mob, opt, OP)
   nc = G.cells.num;
   nf = G.faces.num;

   [press, fp, fp] = impesAssembleStateVars(state, opt.bc, ...
                                            opt.wells, opt.wdp);       %#ok

   N = tpfaExtendedConnections(G, opt.bc, opt.wells);
   i = N ~= 0;
   j = find(~i);

   p    = zeros(size(N));
   p(i) = press(N(i));
   p(j) = fp(mod(j - 1, size(N,1)) + 1);

   dp = p(:,1) - p(:,2);
   g  = sparse([OP.connno ; OP.wconnno], ...
               [OP.cellno ; OP.wcellno], ...
               [OP.gpot   ; OP.wgpot  ] .* ...
               [OP.sgn    ; OP.wsgn]  ) * rho(1 : nc, :);

   dh = repmat(dp, [1, size(g, 2)]) + g;

   flux = trans .* sum(mob .* dh, 2);

   state.flux = flux(1 : nf);

   if ~isempty(opt.wells),
      nw  = numel(opt.wells);
      off = nf;

      for w = 1 : nw,
         np = numel(opt.wells(w).cells);

         assert (numel(state.wellSol(w).flux) == np, ...
                 'Well state is not consistent with well definition');

         state.wellSol(w).flux = flux(off + (1 : np));

         off = off + np;
      end
   end
end
%--------------------------------------------------------------------------

function [state, ok] = update_masses_flux(z0, flux, Af, conmob, dt, pv, OP);%update_masses_flux(state, dp, fflux_dp, gflux, dt, pv, OP)
   %comp_flux = bsxfun(@times, fflux_dp, dp) + gflux;
   %face_flux = bsxfun(@times,state.flux,fmob);
   %well_flux = bsxfun(@times,state.wellSol.flux,wmob);
   %conn_flux = [face_flux(OP.connno);well_flux(OP.wconno)];
   %comp_conflux=Af*conn_flux;
   comp_conflux=mmultiply(Af,bsxfun(@times,flux,bsxfun(@rdivide,conmob,sum(conmob,2))));
   comp_conflux(OP.connno(~OP.active),:)=0;
   %comp_flux =
   %comp_flux(OP.connno(~OP.active),:) = 0;

   cellno = [OP.cellno; OP.wcellno];

   ccontrib = sparse(cellno, 1:numel(cellno), [OP.sgn; OP.wsgn]) ...
              *                                                  ...
              comp_conflux([OP.connno; OP.wconnno], :);
              %comp_flux([OP.connno; OP.wconnno], :);

   dz = - bsxfun(@rdivide, ccontrib, pv);

   state.z = z0 + dt.*dz;

   neg_mass = sum(state.z < 0);
   if any(neg_mass > 0),
      ok = false;
      neg_dist = sprintf(' %d', neg_mass);
      warning('Mass update produces negative masses [%s ] cells', neg_dist);
   else
      ok=true;
   end
end
%--------------------------------------------------------------------------

function [state, dt] = update_masses(state, dp, fflux_dp, gflux, dt, pv, OP)
   comp_flux = bsxfun(@times, fflux_dp, dp) + gflux;
   comp_flux(OP.connno(~OP.active),:) = 0;

   cellno = [OP.cellno; OP.wcellno];

   ccontrib = sparse(cellno, 1:numel(cellno), [OP.sgn; OP.wsgn]) ...
              *                                                  ...
              comp_flux([OP.connno; OP.wconnno], :);

   dz = - bsxfun(@rdivide, ccontrib, pv);

   state.z = state.z + dt.*dz;

   neg_mass = sum(state.z < 0);
   if any(neg_mass > 0),
      neg_dist = sprintf(' %d', neg_mass);
      error('Mass update produces negative masses [%s ] cells', neg_dist);
   end
end

%--------------------------------------------------------------------------

function [A, dA, varargout] = fluid_matrix(p, z, fluid)
   %c    rho  mu  u  R  B  A  dA
   [rho, rho, A,  A, A, A, A, dA] = fluid.pvt(p, z);  %#ok

   if nargout > 2,
      varargout{1} = rho;
   end
end

%--------------------------------------------------------------------------

function F = ...
      residual_cell(luAc, state0, G, fcontrib, dt, pv, opt, OP, luAc_tmp)

   cellno = [OP.cellno(OP.active); OP.wcellno];
   connno = [OP.connno(OP.active); OP.wconnno];
   sgn    = [OP.sgn(   OP.active); OP.wsgn   ];

   if ~isempty(opt.bc),
      f = opt.bc.face;
      c = double(sum(G.faces.neighbors(f,:), 2));
      i = OP.hf(f, c);

      assert (~any(OP.active(i)), ...
              'Internal error defining active cell connections.');

      cellno = [cellno; OP.cellno(i)];
      connno = [connno; OP.connno(i)];
      sgn    = [sgn   ; OP.sgn(   i)];
   end

   ccontrib = sparse(cellno, 1:numel(cellno), sgn) ...
              *                                    ...
              fcontrib(connno, :);

   F = pv .* (1 - sum(solve_single_sys(luAc, state0.z), 2));
   F = F + dt.*sum(solve_single_sys(luAc_tmp, ccontrib), 2);
end

%--------------------------------------------------------------------------

function [F, fac] = residual_well(W, wdp, state, sdyntrans, rdyntrans)
   is_bhp = strcmpi('bhp', { W.type });

   if all(is_bhp),
      % Nothing to do.  All wells controlled by BHP target.
      F = [];
      fac=[];
   else
      nperf = reshape(cellfun('prodofsize', { W.cells }), [], 1);
      pick  = reshape(~is_bhp, [], 1);

      is_res = strcmpi('resv', { W(pick).type });

      rows  = rldecode(pick, nperf);
      nperf = nperf(pick);

      w  = rldecode(1 : numel(nperf), nperf, 2) .';
      pw = rldecode(vertcat(state.wellSol(pick).pressure), nperf);
      pw = pw + wdp(rows);
      pc = state.pressure(vertcat(W(pick).cells));

      dyntrans = [rdyntrans(rows,:), sdyntrans(rows,:)];

      active_p = abs(vertcat(W(pick).compi)) .' > 0;
      active_p(:, is_res) = true;

      np = reshape(sum(active_p), [], 1);

      phase_id = repmat((1 : size(active_p, 1)).', 1, size(active_p, 2));
      phase_id(:, ~is_res) = phase_id(:, ~is_res) + size(active_p, 1);
      phase_id = phase_id(active_p);

      rpos = cumsum([1; reshape(nperf, [], 1)]);
      cpos = cumsum([1; np]);
      r_ix = rldecode([rpos(1 : end - 1), rpos(2 : end) - 1], np   );
      c_ix = rldecode([cpos(1 : end - 1), cpos(2 : end) - 1], nperf);

      r = reshape(mcolon(r_ix(:,1), r_ix(:,2)), [], 1);
      c = reshape(mcolon(c_ix(:,1), c_ix(:,2)), [], 1);

      vals = dyntrans(sub2ind(size(dyntrans), r, phase_id(c)));
      fac  = accumarray(r, vals);

      F = accumarray(w, fac .* (pw - pc)) - vertcat(W(pick).val);
   end
end

%--------------------------------------------------------------------------

function J = approximate_jacobian(luAc, dAc, G, state, state0, fcontrib, ...
                                  w2c, dt, pv, opt, OP, luAc_tmp, ...
                                  fmob, fluid)

   nc = G.cells.num;
   nw = 0;

   cellno = [OP.cellno ; OP.wcellno];   connno = [OP.connno ; OP.wconnno];
   other  = [OP.other  ; OP.wother ];   sgn    = [OP.sgn    ; OP.wsgn   ];

   % Expand global face contributions Aij Tij [ \lambda_\alpha ]_\alpha
   % to one-sided contacts
   %
   fcontrib = fcontrib(connno, :);

   % Compute flux term:
   %
   %  \Delta t sum_\alpha Ai^{-1}\sum_j Aij Tij [ \lambda_\alpha ]_\alpha
   %
   pcontrib = solve_multiple_sys(luAc_tmp, fcontrib, OP);
   c2w      = sum(pcontrib(numel(OP.active) + 1 : end, :), 2);
   ccontrib = dt .* sum(pcontrib(1 : numel(OP.active), :), 2);

   clear pcontrib

   [i, j, v, hf] = deal([]);

   if ~isempty(opt.bc),
      % Contributions from active boundary connections (those for which an
      % explicit boundary condition is specified).
      %
      c  = double(sum(G.faces.neighbors(opt.bc.face, :), 2));
      hf = OP.hf(opt.bc.face, c);
      d  = ccontrib(hf);

      i = [ i ; c ];   j = [ j ; c ];  v = [ v ; d ];
   end

   active_wells = [];
   if ~isempty(opt.wells),
      is_bhp = strcmpi('bhp' , { opt.wells.type });
      pickw  = ~is_bhp;
      nw     = sum(pickw);

      %         c <-> c       (well contributions)
      i = [ i ; OP.wcellno ];
      j = [ j ; OP.wcellno ];
      v = [ v ; dt .* c2w ];

      if nw > 0,
         pickp   = pickw(OP.wother - G.cells.num);
         wcellno = OP.wcellno(pickp);  wother = OP.wother(pickp);

         %         c -> w            w -> c    w <-> w
         i = [ i ; wcellno         ; wother  ; wother ];
         j = [ j ; wother          ; wcellno ; wother ];
         %v = [ v ; -dt.*c2w(pickp) ; dt*-w2c    ; dt*w2c    ];
         v = [ v ; -dt.*c2w(pickp) ; -w2c    ; w2c    ];

         if  ~all(strcmpi('resv', { opt.wells(pickw).type })),
            pmob         = fmob((G.faces.num + 1) : end, :);
            [iw, jw, vw] = differentiate_well_eq(opt, state, pmob, fluid);

            i = [i ; iw];
            j = [j ; jw];
            v = [v ; dt.*vw];
         end
      end
      active_wells = reshape(~is_bhp, [], 1);
   end

   % Contributions from cell-based dynamics ((d/dp) Ai^{-1}) \sum_j (\dots)
   %
   z = bsxfun(@times, pv, state0.z);

   if ~opt.URC_p,
      press = impesAssembleStateVars(state, opt.bc, opt.wells, opt.wdp);
      dp    = pressure_diff(G, press, opt);
      dp    = dt .* sgn .* dp(connno);

      active = [OP.active; true(size(OP.wsgn))];
      if ~isempty(hf), active(hf) = true; end

      dp    = dp .* double(active);

      % add diagonal term assosiated with transport term
      z     = z - sparse(cellno, 1 : numel(cellno), dp)*fcontrib;
   end

   % Compute compressibility-like term:
   %   sum_\alpha Ac^{-1} (dAc/dp) Ac^{-1} z
   %
   dcontrib = solve_single_sys(luAc, z       );  % d <- Ac\z
   dcontrib = mmultiply       (dAc , dcontrib);  % d <- (dAdp)*d
   dcontrib = solve_single_sys(luAc, dcontrib);  % d <- Ac\d
   dcontrib = sum(dcontrib, 2);

   cellno   = cellno  (OP.active);    other = other(OP.active);
   ccontrib = ccontrib(OP.active);

   % Bulk contributions (active/internal connections only):
   %
   %         N(:,1) <-> N(:,2) ; N(:,1) <-> N(:,1)  ; Compr
   %
   i = [ i ; cellno            ; cellno             ; (1 : nc).' ];
   j = [ j ; other             ; cellno             ; (1 : nc).' ];
   v = [ v ; -ccontrib         ; ccontrib           ; dcontrib   ];

   active_dof = cumsum(double([true([nc, 1]); active_wells]));
   i = active_dof(i);
   j = active_dof(j);

   % Assemble final Jacobian matrix.
   J = sparse(i, j, v, nc + nw, nc + nw);
end

%--------------------------------------------------------------------------

function [i, j, v] = differentiate_well_eq(opt, state, pmob, fluid)
   is_bhp  = reshape(strcmpi('bhp' , { opt.wells.type }), [], 1);
   is_resv = reshape(strcmpi('resv', { opt.wells.type }), [], 1);
   pickw   = ~(is_bhp | is_resv);

   % DISREGARD GRAVITY FOR THE TIME BEING!
   pw = [ state.wellSol(pickw).pressure ] .';
   cw = { opt.wells(pickw).cells };
   nperf = cellfun('prodofsize', cw);

   dp = rldecode(pw, nperf) - state.pressure( vertcat( cw{:} ));
   zw = rldecode(vertcat(opt.wells(pickw).compi), nperf);

   cw = vertcat(cw{:});
   zw(dp<0,:) = state.z(cw(dp<0), :);

   [dA, dA] = fluid_matrix(rldecode(pw, nperf), zw, fluid);

   pmob = pmob(rldecode(pickw, cellfun('prodofsize', { opt.wells.cells })), :);
   dflux = mmultiply(dA, bsxfun(@times, vertcat(opt.wells(pickw).WI) .* dp, pmob));

   i = numel(state.pressure) + rldecode(find(pickw(~is_bhp)), nperf(~is_bhp));
   j = cw;
   compi = vertcat(opt.wells(pickw & ~is_bhp).compi);
   compi = rldecode(compi, nperf(~is_bhp));
   v = sum(dflux .* compi, 2);
end

%--------------------------------------------------------------------------

function grav_term = grav_cap(Af, dyntrans, OP, rho)
   connno    = [OP.connno ; OP.wconnno];
   cellno    = [OP.cellno ; OP.wcellno];
   gpot      = [OP.gpot   ; OP.wgpot  ] .* [OP.sgn ; OP.wsgn];

   nconn     = size(Af , 1) / size(rho, 2);
   ncell     = size(rho, 1);

   grav_term = sparse(connno, cellno, gpot, nconn, ncell) * rho;
   grav_term = mmultiply(Af, bsxfun(@times, dyntrans, grav_term));
end

%--------------------------------------------------------------------------

function dp = pressure_diff(G, press, opt)
   N = tpfaExtendedConnections(G, opt.bc, opt.wells);
   i = all(N ~= 0, 2);

   dp = zeros([size(N, 1), 1]);

   dp(i) = press(N(i,1)) - press(N(i,2));  % Pos flux: N(:,1)->N(:,2)
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
