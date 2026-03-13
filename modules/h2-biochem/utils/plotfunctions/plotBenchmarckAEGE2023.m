function plotBenchmarckAEGE2023(nT,pressureNobact,pressurebact,H2loss,...
    totMassH2Nobact,totMassH2bact,G,statesBact,nbact0)
%
% plotEosPRSW -- Plot the pressure, H2 loss due to microbial activity,
% the total mass of H2 along with time as well as the evolution of the
% microbial population.
%
% SYNOPSIS:
%   plotBenchmarckAEGE2023(nT,pressureNobact,pressurebact,H2loss,...
%   totMassH2Nobact,totMassH2bact,G,statesBact,nbact0)
%
% DESCRIPTION:
%   Generates four figures showing:
%       1) the pressure without and with microbial activity (in bar) during
%       time
%       2) the hydrogen loss due to microbial activity (in %) during time
%       3) the total mass of hydrogen (in kg) without and with microbial
%       activity during time
%       4) the microbial population during time
%
%
% PARAMETERS:
%   pres - Vector of pressure in bar for the x-axis.
%   xliqExp - gaz (H2 or CO2) solubility from experiment.
%   xliqsw - gaz (H2 or CO2) solubility from the SW EoS model.
%   xliqpr - gaz (H2 or CO2) solubility from the PR EoS model.
%   nT - number of time iterations
%   pressureNobact - pressure without microbial activity
%   pressurebact - pressure with microbial activity
%   H2loss - Hydrogen loss due to microbial dynamics
%   totMassH2Nobact - total mass of hydrogen, no microorganisms
%   totMassH2bact - total mass of hydrogen with microbial activity
%   G - Grid and Rock Properties
%   statesBact - states during time.
%   nbact0 - Initial microbial concentration (scalar, dimensionless/normalized).
% RETURNS:
%   none (generates figure).
%
% SEE ALSO:
%





% Author: [Stéphanie Delage Santacreu]
% Date: [16/09/2025]
% Organization: [Université de Pau et des Pays de l'Adour, E2S UPPA, CNRS, LFCR, UMR5150, Pau, France]
% ---------------------------------------------------------------------------

% Pressure comparison plot
f01 = figure('Name', 'PressureComparBactMsalt0', 'NumberTitle', 'off');
f01.Position(3:4) = [900, 700];

plot(1:nT, pressurebact./1e5, 'b', 'LineWidth', 2);
hold on;
plot(1:nT, pressureNobact./1e5, 'r--', 'LineWidth', 2);

title('Mean pressure in the pure water reservoir', 'FontSize', 16, 'FontWeight', 'bold');
xlabel('Time (days)', 'FontWeight', 'bold');
ylabel('Pressure (bar)', 'FontWeight', 'bold');

legend({'Pressure, with archae', 'Pressure, no archae'}, ...
    'FontSize', 16, 'Location', 'best');

ax = gca;
ax.FontSize = 16;

%H2 loss
f20=figure('Name','H2loss','NumberTitle','off');
f20.Position(3:4) = [900 700];
plot(1:nT,H2loss,'b','MarkerSize',7,'LineWidth',2)
title(' H2 loss','FontSize',16,'FontWeight','bold','Color','k')
xlabel({'time (days)'},'FontWeight','bold','Color','k')
ylabel({'H2 loss (%)'},'FontWeight','bold','Color','k')
ax = gca;
ax.XMinorTick='on';
ax.FontSize = 16;
legend({'H2 , msalt=0'},...
    'FontSize',16,'TextColor','black','Location','west')

%Total mass of H2
f20=figure('Name','H2totalMass','NumberTitle','off');
f20.Position(3:4) = [900 700];
plot(1:nT,totMassH2Nobact,'b','MarkerSize',7,'LineWidth',2)
hold on;
plot(1:nT,totMassH2bact,'r--','MarkerSize',7,'LineWidth',2)
title(' H2 total mass over time, no salt','FontSize',16,'FontWeight','bold','Color','k')
xlabel({'time (days)'},'FontWeight','bold','Color','k')
ylabel({'H2 mass (kg)'},'FontWeight','bold','Color','k')
ax = gca;
ax.FontSize = 16;
legend({'H2 total mass, no archae','H2 total mass, archae'},...
    'FontSize',16,'TextColor','black','Location','west')


%Microbial population
nbacteria= zeros(nT,1);
ncells=G.cells.num;
for i = 1:nT
    nbacteria(i)=sum(statesBact{i}.nbact);
end

f31=figure('Name','nbacteria','NumberTitle','off');
f31.Position(3:4) = [900 700];
plot(0:nT,[ncells*nbact0;nbacteria],'b','MarkerSize',7,'LineWidth',2)
title('Total methanogenic Archae population','FontSize',16,'FontWeight','bold','Color','k')
xlabel({'time (days)'},'FontWeight','bold','Color','k')
ylabel({'N_{archae}'},'FontWeight','bold','Color','k')
ax = gca;
ax.FontSize = 16;
legend({'N_{archae}, msalt=0'},'FontSize',16,'TextColor','black',...
    'Location','best')


end


%{
Copyright 2009-2025 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify it under the terms of the GNU 
General Public License as published by the Free Software Foundation, either version 3 of 
the License, or (at your option) any later version.

MRST is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without 
even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with MRST. If not, 
see <http://www.gnu.org/licenses/>.
%}