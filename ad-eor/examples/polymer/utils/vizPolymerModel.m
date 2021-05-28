%% Script to visualize some of the surfactant models properties
%

G     = model.G;
rock  = model.rock;
fluid = model.fluid;


%% Plotting the grid, wells, and porosity and permeability

switch example_name
  case 'blackoil2D'
    h=figure(); clf
    set(h, 'Position', [100, 100, 900, 600]);
    subplot(2,2,1);
    plotGrid(G, 'FaceColor', 'none');
    % plot wells, blue for injector, red for producer.
    W = schedule.control(1).W; sgn = [W.sign];
    plotWell(G, W(sgn>0), 'color', 'b')
    plotWell(G, W(sgn<0), 'color', 'r')
    axis tight off; view(20,8);

    prop = {'PERMX', 'PERMZ', 'porosity'};
    for i = 1:3
        if strcmp(prop{i},'PERMX')
            s_plot = rock.perm(:,1)/(milli*darcy);
        end
        if strcmp(prop{i},'PERMZ')
            s_plot = rock.perm(:,3)/(milli*darcy);
        end
        if strcmp(prop{i},'porosity')
            s_plot = rock.poro;
        end
        subplot(2, 2, i+1);
        plotCellData(G, s_plot, 'EdgeColor','k');
        title(prop{i});
        axis tight off; view(20,8); colorbar;
    end
  case 'spe10'
    figure; clf;
    W = schedule.control(1).W;
    sgn = [W.sign];
    % plot injector wells
    plotWell(G, W(sgn>0), 'fontsize', 0, 'color', 'b');
    % plot producer wells
    plotWell(G, W(sgn<0), 'fontsize', 0, 'color', 'r');
    plotCellData(G, log10(rock.perm(:, 1)),'EdgeColor','none')
    logColorbar();
    title('Permeability (x-direction), spe10')
    set(gca,'dataasp',[60 60 1]);
    axis tight,
    view(80, 80);
end


%% Plotting the inital saturations and pressures
switch example_name
  case 'blackoil2D'
    h=figure(); clf
    set(h, 'Position', [100, 100, 900, 600]);
    % plotting the initial saturations
    phases = {'water','oil','gas'};
    for i=1:3
        subplot(2,2,i+1);
        plotCellData(G, state0.s(:,i), 'EdgeColor','k');
        title(['Initial ', phases{i}, ' saturation'])
        axis tight off; view(20,8); caxis([0, 1]); colorbar;
    end
    % plotting the inital pressure
    subplot(2, 2, 1);
    plotCellData(G, state0.pressure/barsa, 'EdgeColor','k');
    title('Inital pressure (Bar)');
    axis tight off; view(20,8); colorbar;
end

%% Plotting the polymer properties
h=figure(); clf
set(h, 'Position', [100, 100, 1200, 300]);

% plotting the viscosity multiplier, which describes the viscosity
% enhancement effect due to polymer mixing
subplot(1,3,1);
c = (0:0.3:3);
vismult = fluid.muWMult(c);
plot(c, vismult, '-*', 'linewidth', 2);
xlabel('polymer concentration (kg/m^3)');
ylabel('viscosity multiplier');
title('Polymer viscosity multiplier');

% plotting the polymer adsportion function
subplot(1,3,2);
ads = fluid.ads(c);
plot(c, ads, '-*', 'linewidth', 2);
xlabel('polymer concentration (kg/m^3)');
ylabel('adsorped polymer (kg/kg)');
title('Polymer adsorption');

% plotting the shear factor
% from the plot we can see, the polymer in this example holds
% shear-thinning flow rheology.
if model.usingShear
    subplot(1,3,3);
    v = 10.^(-11:1:0);
    shear_factor = fluid.plyshearMult(v);
    semilogx(v, shear_factor, '-*', 'linewidth', 2);
    xlabel('water velocity (m/s)');
    ylabel('shear factor');
    title('Shear effect');
end


drawnow;

% <html>
% <p><font size="-1">
% Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.
% </font></p>
% <p><font size="-1">
% This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
% </font></p>
% <p><font size="-1">
% MRST is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% </font></p>
% <p><font size="-1">
% MRST is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% </font></p>
% <p><font size="-1">
% You should have received a copy of the GNU General Public License
% along with MRST.  If not, see
% <a href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses</a>.
% </font></p>
% </html>
