%% Script to visualize some of the surfactant models properties
%

G     = model.G;
rock  = model.rock;
fluid = model.fluid;


%% Grid plot with wells
%

figure
clf
W = schedule.control(1).W;
sgn = [W.sign];
% plot injector wells
plotWell(G, W(sgn>0), 'fontsize', 0, 'color', 'b');
% plot producer wells
plotWell(G, W(sgn<0), 'fontsize', 0, 'color', 'r');
switch example_name
  case '1D'
    plotGrid(G, 'FaceColor', 'none')
    view(3), axis equal
  case '2D'
    plotGrid(G, 'FaceColor', 'none')
    view(70, 30);
  case 'spe9'
    plotCellData(G, log10(rock.perm(:, 1)),'EdgeColor','none')
    logColorbar();
    title('Permeability (x-direction), spe10')
    set(gca,'dataasp',[60 60 1]);
    axis tight,
    view(30, 86);
  case 'spe10'
    plotCellData(G, log10(rock.perm(:, 1)),'EdgeColor','none')
    logColorbar();
    title('Permeability (x-direction), spe10')
    set(gca,'dataasp',[60 60 1]);
    axis tight,
    view(80, 80);
end


%% Relative permeabilities
% The relative permeabilities depend on the saturation and surfactant
% concentration. The dependence with respect to the surfactant concentration is
% computed by interpolation between two relative permeability curves, one
% without surfactant and one with a fully mixed saturated solution with
% surfactant. See the document ad-eor/docs/surfactant_model.pdf for more detail.
%
% The interpolation factor depends on the capillary number (N_c), which measures the
% relative importance of the interfacial tension with respect to the viscous
% forces. See the functions computeCapillaryNumber.m, computeRelPermSft.m. As
% an example, we choose N_c =
%
% From those curves, we can see in particular how the surfactant mobilize the
% residual oil.
%

s = (0 : 0.05 : 1)';

switch example_name
  case {'1D', '2D'}
    % We use deck input. For such input, different relative permeabilities can be
    % assigned to different cell. In our case, they are the same everywhere and
    % we choose the first cell to comput them for plotting.
    krW = fluid.krW(s, 'cellInx', 1);
    krOW = fluid.krOW(s, 'cellInx', 1);
    krWSft = fluid.krWSft(s, 'cellInx', 1);
    krOWSft = fluid.krOWSft(s, 'cellInx', 1);
    if all(isfield(fluid,{'krG','krOG'}))
       krG = fluid.krG(s, 'cellInx', 1);
       krOG = fluid.krOG(s, 'cellInx', 1);
    end
  case 'spe9'
    krW = fluid.krW(s, 'cellInx', 1);
    krOW = fluid.krOW(s, 'cellInx', 1);
    krWSft = fluid.krWSft(s, 'cellInx', 1);
    krOWSft = fluid.krOWSft(s, 'cellInx', 1);
    if all(isfield(fluid,{'krG','krOG'}))
       krG = fluid.krG(s, 'cellInx', 1);
       krOG = fluid.krOG(s, 'cellInx', 1);
    end
  case 'spe10'
    krW = fluid.krW(s);
    krOW = fluid.krOW(s);
    krWSft = fluid.krWSft(s);
    krOWSft = fluid.krOWSft(s);
end

if all(isfield(fluid,{'krG','krOG'}))
    
    figure()
    plot(s, krW, s, krWSft, 'LineWidth',2);
    xlabel('Water saturation');
    ylabel('Relative permeability');
    title('Water Relative Permeabilities');
    legend({'$kr_w^{\textrm{nosurf}}$ (no surfactant)', ['$kr_w^{\textrm{surf}}$ ' ...
        '(fully-mixed with surfactant)']}, 'interpreter', 'latex', ...
        'location', 'best')
    
    figure()
    plot(1 - s, krOW, 1 - s, krOWSft,'LineWidth',2);
    xlabel('Water saturation');
    ylabel('Relative permeability');
    title('Oil Relative Permeabilities');
    legend({'$kr_ow^{\textrm{nosurf}}$ (no surfactant)', ['$kr_ow^{\textrm{surf}}$ ' ...
        '(fully-mixed with surfactant)']}, 'interpreter', 'latex', ...
        'location', 'best')
    
    figure()
    plot(s, krG, 1 - s, krOG,'LineWidth',2);
    xlabel('Gas saturation');
    ylabel('Relative permeability');
    title('Gas Relative Permeabilities');
    legend({'$kr_g$ ', '$kr_o$ '}, 'interpreter', 'latex', ...
        'location', 'best')
else
figure()
plot(s, krW, s, krWSft, 'LineWidth',2);
xlabel('Water saturation');
ylabel('Relative permeability');
title('Water Relative Permeabilities');
legend({'$kr_w^{\textrm{nosurf}}$ (no surfactant)', ['$kr_w^{\textrm{surf}}$ ' ...
                    '(fully-mixed with surfactant)']}, 'interpreter', 'latex', ...
        'location', 'best')

figure()
plot(1 - s, krOW, 1 - s, krOWSft,'LineWidth',2);
xlabel('Water saturation');
ylabel('Relative permeability');
title('Oil Relative Permeabilities');
legend({'$kr_o^{\textrm{nosurf}}$ (no surfactant)', ['$kr_o^{\textrm{surf}}$ ' ...
                    '(fully-mixed with surfactant)']}, 'interpreter', 'latex', ...
       'location', 'best')

end
%% Interfacial tension
% The interfacial tension decreases as the surfactant concentration increases
% data from SPE paper 145036

figure()
c = (0 : 0.1 : 10)';
plot(c, fluid.ift(c),'LineWidth',2);
xlabel('Surfactant concentration');
ylabel('Relative permeability');
title('Interfacial surface tension');

%% Viscosity
% The surfactant, as polymers, modifies the viscosity of the water and makes
% it thicker.
%

figure()
c = (0 : 10 : 200)';
plot(c, convertTo(fluid.muWSft(c), centi*poise), 'LineWidth', 2);
xlabel('Surfactant concentration');
ylabel('Viscosity (cP)');
title('Water Viscosity at reference pressure');

drawnow;

%% Copyright notice

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
