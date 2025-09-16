function reports = postprocessStates(Gt, states, rock, fluid, schedule, traps, dh)
%
% This function takes the results of a simulation ('states', a cell array of
% states for each timestep), and returns a corresponding structure array of
% 'reports'.  The prepared report for a timestep contains the timestep state
% itself, but also other information such as the computed plume heights, and
% the trapping state of the CO2.  
%
% Reports are needed as input to generate inventory plots using
% 'plotTrappingDistribution'.
% 
% Currently, only rate controlled wells are supported (not pressure-controlled).
% 
% SYNOPSIS:
%   function reports = postprocessStates(Gt, states, rock, fluid, schedule, traps, dh)
%
% DESCRIPTION:
%
% PARAMETERS:
%   Gt       - top surface grid used in the simulation
%   states   - result from a simulation (cell array of states, including
%              initial state)
%   rock     - rock object used in the simulation
%   fluid    - fluid object used in the simulation
%   schedule - schedule used in the simulation (NB: only rate controlled
%              wells supported at present)
%   traps    - trapping structure (from trapAnalysis of Gt)
%   dh       - subscale trapping capacity (empty, or one value per grid cell
%              of Gt)
%
% RETURNS:
%   reports - a structure array of 'reports', one per timestep in the 
%             simulation.
%
% SEE ALSO:
%   `plotTrappingDistribution`.

%{
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

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
   
   assert( numel(states) == numel(schedule.step.val)+1 , ...
       'Ensure the initial state has been included in the varargin ''states''.')

   for i = 1:numel(states)

      if i == 1
         % initial state can contain non-zero co2 saturations (or heights)
         ntg = ones(Gt.cells.num,1);
         if isfield(rock,'ntg')
            ntg = rock.ntg;
         end
         tot_inj = Gt.cells.volumes .* Gt.cells.H .* rock.poro .* ntg .* (1-fluid.res_water) ...
             .* states{1}.s(:,2) .* fluid.rhoG(states{1}.pressure);
         % tot_inj = Gt.cells.volumes .* rock.poro .* ntg .* (1-fluid.res_water) ...
         %     .* h .* fluid.rhoG(states{1}.pressure);
         tot_inj = sum(tot_inj);
         reports(i).t = 0; %#ok
         reports(i).W = []; %#ok
      else
         reports(i).t = sum(schedule.step.val(1:i-1)); %#ok
         reports(i).W = schedule.control(schedule.step.control(i-1)).W; %#ok
         
         assert(all(cellfun(@(x) strcmpi(x, 'rate'), {reports(i).W.type})));
         tot_inj = tot_inj + sum([reports(i).W.val]) * schedule.step.val(i-1) * fluid.rhoGS;
      end
      
      reports(i).sol       = states{i}; %#ok
      
      rs = 0;
      if isfield(states{i}, 'rs')
         rs = states{i}.rs;
      end
      if isfield(states{i}, 'sGmax')
          smax = states{i}.sGmax; % we operate with dissolution
      else
          smax = states{i}.smax(:,2); % no dissolution.  Current max = historical max
      end
      reports(i).masses = massTrappingDistributionVE(states{i}.s(:,2), smax, states{i}.pressure, rs, Gt, fluid, rock, traps, dh);
      
      leaked = tot_inj - sum(reports(i).masses);
      reports(i).tot_inj = tot_inj;
      reports(i).masses = [reports(i).masses, leaked]; %#ok
   end
   
end
