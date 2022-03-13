function compareWellrates(initSchedule, optimisedSchedule, co2RefRho)
%Undocumented Utility Function

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

   init_rates = [initSchedule.control(1).W.val];
   opt_rates  = [optimisedSchedule.control(1).W.val];

   val = [init_rates; opt_rates]';
    
    % values are currently in m3/sec using the CO2 reference density.  We
    % convert it to Mt/year
    val = val * co2RefRho * year / 1e6 / 1e3;
    
    figure; bar(val);
    set(gca, 'xlim', [0 size(val, 1) + 1]);
    xlabel('Well', 'fontsize', 14);
    ylabel('Rate (Mt/year)', 'fontsize', 14);
    set(gca, 'fontsize', 14); % increase fontsize on axes.
    box off;
end
