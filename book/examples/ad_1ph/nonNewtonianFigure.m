mu0 = 100*centi*poise;
for nsim=1:4
   switch nsim
      case 1
         fluidModel = struct('mu0', mu0, 'nmu',  1, 'Kc', .1);
         nonNewtonianCell;
      case 2
         fluidModel = struct('mu0', mu0, 'nmu', .3, 'Kc', .1);
         nonNewtonianCell;
      case 3
         fluidModel = struct('mu0', mu0, 'nmu', .3, 'Kc', .1);
         wellAvg = true;
         nonNewtonianFace;
      case 4
         fluidModel = struct('mu0', mu0, 'nmu', .3, 'Kc', .1);
         wellAvg = false;
         nonNewtonianFace;
   end
   avgpres(:,nsim) = mean([sol(2:end).pressure]/barsa);        %#ok<SAGROW>
   rate   (:,nsim) = [sol(2:end).qS]*day;                      %#ok<SAGROW>
   time   (:,nsim) = [sol(2:end).time]/day;                    %#ok<SAGROW>
   mineta (:,nsim) = etamin;                                   %#ok<SAGROW>
   minweta(:,nsim) = etawmin;                                  %#ok<SAGROW>
   meaneta(:,nsim) = etamean;                                  %#ok<SAGROW>
end


%%
figure('Position', [440 375 840 420]);
subplot(1,2,1);
stairs(time, rate, 'LineWidth', 2);
set(gca,'FontSize',12);
xlabel('time [days]'); ylabel('rate [m^3/day]');
legend('Newtonian', 'Cell-based', 'Face-based', 'Face-based (not well)', ...
       'Location', 'NorthEast');
subplot(1,2,2);
plot(time, avgpres,'o','LineWidth', 2);
set(gca,'FontSize',12,'YAxisLocation','right');
xlabel('time [days]'); ylabel('avg pressure [bar]');
legend('Newtonian', 'Cell-based', 'Face-based', 'Face-based (not well)', ...
       'Location', 'SouthEast');
set(gcf,'PaperPositionMode','auto');

%%
subplot(1,2,2), cla
plot(time,mineta,'-','LineWidth',2);
hold on
plot(time,minweta,'--','LineWidth',2);
plot(time,meaneta,'-.','LineWidth',2);
hold off
set(gca,'FontSize',12,'YAxisLocation','right');
xlabel('time [days]'); ylabel('shear multiplicator [1]');

%{
Copyright 2009-2018 SINTEF Digital, Mathematics & Cybernetics.

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
