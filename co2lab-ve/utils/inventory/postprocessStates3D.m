function reports = postprocessStates3D(G, states, rock, fluid, schedule, res_water, res_gas)

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
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

    % convert information to VE, and call `postprocessStates`
    Gt = topSurfaceGrid(G);
    rock2D = averageRock(rock, Gt);
    ta = trapAnalysis(Gt, false);
    statesVE = states2VE(states, Gt, fluid, rock.poro);
    
    fluidVE.res_water = res_water;
    fluidVE.res_gas = res_gas;
    fluidVE.rhoGS = fluid.rhoGS;
    fluidVE.rhoWS = fluid.rhoWS;
    fluidVE.rhoG = @(p, varargin) fluid.bG(p, varargin{:}) * fluidVE.rhoGS;
    fluidVE.bG = fluid.bG;
    fluidVE.bW = fluid.bW;

    reports = postprocessStates(Gt, statesVE, rock2D, fluidVE, schedule, ta, []);
    
    % add-in original 3D states
    for i=1:numel(states)
        reports(i).sol3D = states{i};
    end
    
end
