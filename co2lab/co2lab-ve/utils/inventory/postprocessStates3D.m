function reports = postprocessStates3D(G, states, rock, fluid, schedule, ...
                                       res_water, res_gas, varargin)

% This function is a wrapper that calls 'postprocessStates' on a result from a 3D
% simulation. It takes the result of a simulation ('states', a cell array of
% states for each timestep), and returns a corresponding structure array of
% 'reports'.  The prepared report contains each timestep state, but also
% other information such as the computed plume heights and the trapping state
% of the CO2.
%
% Reports are needed as input to generate inventory plots using
% 'plotTrappingDistribution'. 
% 
% Currently, only rate controlled wells are supported (not pressure-controlled).
% 
% SYNOPSIS:
%   function reports = postprocessStates3D(G, states, rock, fluid, schedule, traps,
%   res_water, res_gas)
%
% DESCRIPTION:
%
% PARAMETERS:
%   G         - 3D simulation grid used in the simulation
%   states    - result from the 3D simulation (cell array of states, including
%               initial state)
%   rock      - rock object used in the simulation
%   fluid     - fluid object used in the simulation
%   schedule  - schedule used in the simulation (NB: only rate controlled
%               wells supported at present)
%   res_water - residual water saturation
%   res_gas   - residual co2 saturation
%
% OPTIONAL PARAMETER
%   traps     - If a trap analysis has already been carried out on the top 
%               surface grid, it can be directly provided here (otherwise 
%               it will be computed internally).
% RETURNS:
%   reports - a structure array of 'reports', one per timestep in the 
%             simulation.
%
% SEE ALSO:
%   `plotTrappingDistribution`, `postprocessStates`

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
GNU General Public License for more details.b

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
    opt.traps = [];
    opt = merge_options(opt, varargin{:});

    for i = 1:numel(states)
    
        if i == 1
            % initial state can contain non-zero co2 saturations (or heights)
            ntg = ones(G.cells.num,1);
            if isfield(rock,'ntg')
                ntg = rock.ntg;
            end
            tot_inj = G.cells.volumes .* rock.poro .* ntg .* states{1}.s(:,2) ...
                      .* fluid.bG(states{1}.pressure) .* fluid.rhoGS;
            tot_inj = sum(tot_inj);
            reports(i).t = 0; %#ok
            reports(i).W = []; %#ok
        else
            reports(i).t = sum(schedule.step.val(1:i-1)); %#ok
            reports(i).W = schedule.control(schedule.step.control(i-1)).W; %#ok
            
            assert(all(cellfun(@(x) strcmpi(x, 'rate'), {reports(i).W.type})));
            tot_inj = tot_inj + sum([reports(i).W.val]) * schedule.step.val(i-1) * fluid.rhoGS;
            
        end
        
        reports(i).sol3D = states{i}; %#ok

        rs = 0;
        if isfield(states{i}, 'rs')
            rs = states{i}.rs;
        end

        smax = states{i}.s(:,2);
        if isfield(states{i}, 'sMax')
            smax = states{i}.sMax(:, 2);
        end
        
        reports(i).masses = ...
            massTrappingDistribution3D(states{i}.s(:,2), smax, ...
                                       states{i}.pressure, rs, G, fluid, ...
                                       rock, res_water, res_gas, 'trapstruct', opt.traps);
        
        leaked = tot_inj - sum(reports(i).masses);
        reports(i).tot_inj = tot_inj;
        reports(i).masses = [reports(i).masses, leaked]; %#ok
        
    end
    
end
