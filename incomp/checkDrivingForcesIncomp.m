function do_solve = checkDrivingForcesIncomp(g, opt)
%Check if a solution is necessary for incompressible solvers. 
%
% SYNOPSIS:
%   W = addWell(W, G, rock, cellInx)
%   W = addWell(W, G, rock, cellInx, 'pn', pv, ...)
%
% REQUIRED PARAMETERS:
%   g       - Valid MRST grid intended for simulation.
%
%   opt     - Merged options struct with standard incompressible options.
%             Specifically, the forces supported by incompTPFA and
%             incompMimetic should be included.
%
%
% RETURNS:
%
%   do_solve - Logical indicating if forces are well specified and should
%              be solved. This check verifies that the problem is
%              meaningful for incompressible flow.
% SEE ALSO:
%   `incompMimetic`, `incompTPFA`.

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

   pressure_bc = ~isempty(opt.bc) && ...
                 any(strcmpi('pressure', opt.bc.type));
   well_bhp    = ~isempty(opt.wells) && ...
                 any(strcmpi('bhp', { opt.wells.type }));
   if ~(pressure_bc || well_bhp)
      sum_rates = 0;
      norm_rates = 0;
      if ~isempty(opt.wells)
          rates = vertcat(opt.wells.val);
          sum_rates = sum(rates);
          norm_rates = max(norm_rates, norm(rates, inf));
      end
      if ~isempty(opt.src)
          sum_rates = sum_rates + sum(opt.src.rate);
          norm_rates = max(norm_rates, norm(opt.src.rate, inf));
      end
      if ~isempty(opt.bc)
          sum_rates = sum_rates + sum(opt.bc.value);
          norm_rates = max(norm_rates, norm(opt.bc.value, inf));
      end
      if abs(sum_rates) > 1e-12
         warning('incomp:MassBalance:NotFulfilled', ...
              ['Well rates and flux bc must sum up to 0 \n', ...
               'when there are no bhp constrained wells or pressure bc.'...
               'Results may not be reliable.\n']);
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
                            isempty(opt.bcp)
                            ]);
   % Assemble (and solve) system even in absence of external driving forces
   % if the caller requested 'MatrixOutput'.
   do_solve = do_solve || opt.MatrixOutput;

   if isfield(opt, 'rhs')
       do_solve = do_solve || ~isempty(opt.rhs);
   end
   if ~do_solve,
      warning('incomp:DrivingForce:Missing',                   ...
             ['No external driving forces present in model--', ...
              'state remains unchanged.\n']);
   end

end
