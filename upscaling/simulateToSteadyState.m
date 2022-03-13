function [state, report] = simulateToSteadyState(state, G, rock, fluid, dt, varargin)
% Simulate until steady state is reached
%
% SYNOPSIS:
%    [state, report] = simulateToSteadyStatePeriodic(state, G, rock, fluid, dt)
%    [state, report] = simulateToSteadyStatePeriodic(state, G, rock, fluid, dt, 'pn', pv,...)
%
% PARAMETERS:
%   state  - Reservoir and well solution structure either properly
%            initialized from functions 'initResSol' and 'initWellSol'
%            respectively, or the results from a previous call to function
%            'incompTPFA' and, possibly, a transport solver such as
%            function 'implicitTransport'.
%
%   G      - Valid grid structure, see grid_structure.
%
%   rock   - Rock containing fields 'poro' and 'perm'.
%
%   dt     - Initial time step used while trying to reach steady state. The
%            time steps will be increased if the changes are miniscule
%            during each timestep.
%
%   fluid  - Fluid object.
%
%   Contains
%
%   'pn'/pv - List of 'key'/value pairs defining optional parameters.  The
%             supported options are:
%
%   fluid_nc - Fluid without capillary pressure. Will remove capillary
%              pressure from fluid argument if not provided.
%
%   psolver  - @(state,fluid) function handle to pressure solver.
%
%   bc       - Boundary conditions
%
%   bcp      - Periodic boundary conditions (if the grid is periodic)
%
%   solve_pressure - Pressure is solved each iteration
%
%   dt_max   - Maximum time step allowed.
%
%   diff_tol - The tolerance used for checking if stationary state has been
%              found.
%
%   <other>  - Passed onto implicitTransport. Look there for explanations
%              of the various parameters.
% RETURNS:
%   state    - Steady state.
%
% COMMENTS:
%
%
% SEE ALSO:
%

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

   try
      require mimetic
   catch
      mrstModule add mimetic
   end

opt =struct('fluid_nc',    [],...
            'psolver',     [],... % @(state, fluid) format
            'verbose',     mrstVerbose,...
            'bc',          [],...
            'bcp',         [],...
            'solve_pressure', true,...
            'dt_max',      300*year,...
            'diff_tol',    1e-5,...
            'nltol'   ,    1.0e-6, ...  % Non-linear residual tolerance
            'lstrials',    10    , ...  % Max no of line search trials
            'maxnewt' ,    20    , ...  % Max no. of NR iterations
            'tsref'  ,     1     , ...  % Time step refinement
            'max_it',      100,...
            'max_newton',  200,...
            'trans',       [],...
            'dhfz',        []...
             );
   opt = merge_options(opt, varargin{:});

   assert(xor(isempty(opt.bc), isempty(opt.bcp)), 'We want either bc or periodic bc - not both!');
   if ~isempty(opt.bc)
       bc_active = opt.bc;
   else
       bc_active = opt.bcp;
   end
   if isempty(opt.fluid_nc)
       % Strip capillary pressure from fluid object if it exists and
       % another fluid without capillary pressure isn't provided.
       if isfield(fluid, 'pc')
           fluid_nc = rmfield(fluid, 'pc');
       else
           fluid_nc = fluid;
       end
   end

   if isempty(opt.psolver)
       if isempty(opt.bcp)
           S = computeMimeticIP(G, rock);
           opt.psolver = @(state, fluid, dummy) ...
              incompMimetic(state, G, S, fluid, 'bc', opt.bc);
       else
           error(['Please provide pressure solver through option psolver'...
               ' with a @(state, fluid, bc/bcp) format']);
       end
   end
   psolver = @(state) opt.psolver(state, fluid_nc, bc_active);

   state = psolver(state);

   transport_cost=0;
   it = 0;
   ds = 0;

   init_state = state;
   stationary = false;

   while (~stationary && it < opt.max_it && transport_cost < opt.max_newton)
       state_old=state;

       [state,report]=implicitTransport(state, G, dt, rock, fluid,...
                                        'Trans'     ,   opt.trans,...
                                        'verbose'   ,   opt.verbose,...
                                        'nltol'     ,   opt.nltol, ...     % Non-linear residual tolerance
                                        'lstrials'  ,   opt.lstrials, ...  % Max no of line search trials
                                        'maxnewt'   ,   opt.maxnewt, ...   % Max no. of NR iterations
                                        'tsref'     ,   opt.tsref, ...     % Time step refinement
                                        'init_state',   init_state, ...
                                        'dhfz'      ,   opt.dhfz,...
                                        'bc', opt.bc);
      ds_prev = ds;
      ds      = norm(state_old.s-state.s,inf);
      transport_cost = transport_cost+report.iterations;

      dispif(opt.verbose, ['Iteration ',num2str(it),' DT ', num2str(dt/year),' Error',num2str(ds_prev/(1-ds/ds_prev))])
      dispif(opt.verbose,['Newton iteration ',num2str(report.iterations)])
      dispif(opt.verbose,['Iteration ',num2str(it),' DT ', num2str(dt/year),' Error ', num2str(ds)])

      if(report.success)
          % If the step was successful, reset initial state and increase
          % possibly increase timestep
         init_state=state;
         if(opt.solve_pressure)
            state=psolver(state);
         end
         if(norm(state_old.s-state.s,inf)<1e-2 && norm((state_old.flux-state.flux)/max(abs(state.flux)),inf)<1e-1)
            dt=min(opt.dt_max,dt*2);
         end

         if((ds*year/dt)<opt.diff_tol)
            stationary=true;
         end
      else
         dispif(opt.verbose, 'Cutting time step');
         init_state.s=state_old.s;
         dt=dt/2;
      end
      it=it+1;
   end
   if opt.verbose
       disp(['Total newton iterations is ', num2str(transport_cost)])
       disp(['Pressure solves are ', num2str(it)])
       disp(['Largest step was ' num2str(dt/year),' year'])
       disp(['Change per year was ' num2str((ds*year)/dt)])
   end
   report=struct('stationary',stationary,'steps',it,'newton_steps',transport_cost);
end
