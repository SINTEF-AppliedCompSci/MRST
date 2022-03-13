function reports = makeReports(Gt, states, rock, fluid, schedule, residual, traps, dh)
%
% This function does intermediate processing of simulation data in order to
% generate inventory plots using 'plotTrappingDistribution'.
% 
% Currently, only rate controlled wells are supported (not pressure-controlled).
% 
% SYNOPSIS:
%   function reports = makeReports(Gt, states, rock, fluid, schedule, residual, traps, dh)
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
%   residual - residual saturations, on the form [s_water, s_co2]
%   traps    - trapping structure (from trapAnalysis of Gt)
%   dh       - subscale trapping capacity (empty, or one value per grid cell
%              of Gt)
%
% RETURNS:
%   reports - a structure array of 'reports', that can be provided to the
%   'plotTrappingDistribution' function.
%
% SEE ALSO:
%   `plotTrappingDistribution`.

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
   
   assert( numel(states) == numel(schedule.step.val)+1 , ...
       'Ensure the initial state has been included in the varargin ''states''.')
   
   for i = 1:numel(states)

      [h, h_max] = compute_plume_height(Gt, states{i}, residual(1), residual(2));

      if i == 1
         % initial state can contain non-zero co2 saturations (or heights)
         ntg = ones(Gt.cells.num,1);
         if isfield(rock,'ntg')
            ntg = rock.ntg;
         end
         tot_inj = Gt.cells.volumes .* rock.poro .* ntg .* (1-fluid.res_water) ...
             .* h .* fluid.rhoG(states{1}.pressure);
         tot_inj = sum(tot_inj);
         reports(i).t         = 0; %#ok
         reports(i).W         = []; %#ok
      else
         reports(i).t = sum(schedule.step.val(1:i-1)); %#ok
         reports(i).W = schedule.control(schedule.step.control(i-1)).W; %#ok
         
         assert(all(cellfun(@(x) strcmpi(x, 'rate'), {reports(i).W.type})));
         tot_inj = tot_inj + sum([reports(i).W.val]) * schedule.step.val(i-1) * fluid.rhoGS;
      end
      
      reports(i).sol       = states{i}; %#ok
      reports(i).sol.h     = h; %#ok
      reports(i).sol.h_max = h_max; %#ok
      
      rs = 0;
      if isfield(states{i}, 'rs')
         rs = states{i}.rs;
      end
      
      reports(i).masses    = massTrappingDistributionVEADI(Gt                   , ...
                                                        reports(i).sol.pressure , ...
                                                        reports(i).sol.s(:,2)   , ...
                                                        reports(i).sol.s(:,1)   , ...
                                                        reports(i).sol.h        , ...
                                                        reports(i).sol.h_max    , ...
                                                        rock                    , ...
                                                        fluid                   , ...
                                                        traps                   , ...
                                                        dh                      , ...
                                                        'rs', rs); %#ok
      leaked = tot_inj - sum(reports(i).masses);
      reports(i).masses = [reports(i).masses, leaked]; %#ok
   end
   
end
% ----------------------------------------------------------------------------

function [h, h_max] = compute_plume_height(Gt, state, sw, sr)
    
    if isfield(state, 'sGmax')
       smax = state.sGmax; % we operate with dissolution
    else
       smax = state.smax(:,2); % no dissolution.  Current max = historical max
    end
    [h, h_max] = upscaledSat2height(state.s(:,2), smax, Gt, 'resSat', [sw, sr]);
end
