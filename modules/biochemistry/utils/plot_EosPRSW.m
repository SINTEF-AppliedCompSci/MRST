function plot_EosPRSW(pres,xliqExp,xliqsw,xliqpr)
% Function which plot the hydrogen or dioxyde carbone solubility 
% from experiment, SW EoS and PR EoS.

% Author: [Stéphanie Delage Santacreu]
% Date: [16/09/2025]
% Organization: [Université de Pau et des Pays de l'Adour, E2S UPPA, CNRS, LFCR, UMR5150, Pau, France]


patm = 1e5; % Atmospheric pressure in Pa
presbar=pres./patm;
min_xliq=[min(xliqExp),min(xliqsw),min(xliqpr)];
max_xliq=[max(xliqExp),max(xliqsw),max(xliqpr)];

f21=figure('Name','solubility_SwPr','NumberTitle','off');
f21.Position(3:4) = [900 700];

plot(presbar,xliqsw,'k*','MarkerSize',7,'LineWidth',2)
hold on
plot(presbar,xliqpr,'r square','MarkerSize',8,'LineWidth',2)
hold on
plot(presbar,xliqExp,'b o','MarkerSize',8,'LineWidth',2)

title('solubility in pure water PR vs SW','FontSize',16,'FontWeight','normal','Color','k')
xlabel({'pressure (bar)'},'FontWeight','normal','Color','k')
ylabel('molar fraction (-)','FontWeight','normal','Color','k')
ax = gca;
ax.FontSize = 16; 

legend({'SW','PR','Exp'},...
    'FontSize',16,'TextColor','black',...
    'Location','best')
xlim([min(presbar)-10 max(presbar)+10])
ylim([min(min_xliq)-1e-4 max(max_xliq)+1e-4])

end
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