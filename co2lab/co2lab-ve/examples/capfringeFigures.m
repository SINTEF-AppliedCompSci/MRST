% Plot graphs of the following, with and without residual trapping:
% 
% - saturation vs. depth
% - relperm as fct. of upscaled saturation
% - capillary pressure as fct. of upscaled saturation
mrstModule add co2lab-common co2lab-ve

gravity on
%% Make grid

% make a single column grid as a basis for the test; use a very high vertical
% resolution (10 cm thick cells)
G = cartGrid([1, 1, 1000], [1, 1, 50] * meter);
[Gt, G] = topSurfaceGrid(G);

%% Create fluids objects with upscaled saturation, relperm and capillary
%% pressure functions

srw = 0.2; %0.2;
src = 0.25;
rock2D.poro = 0.3;
T = 300;
[C, alpha] = deal(0.4, 0.5); % specifies the capillary pressue function
kr3D = @(s) max((s-src) ./ (1-src), 0) .^ 2;

% fluid with capillary fringe
fluidCP = makeVEFluid(Gt, rock2D, 'P-scaled table', ...
                      'reservoirT', T, ...
                      'residual', [srw, src], ...
                      'invPc3D', [C, alpha], ...
                      'kr3D', kr3D);

% simple, sharp-interface fluid
fluidSI = makeVEFluid(Gt, rock2D, 'sharp_interface_simple', ...
                      'reservoirT', T, ...
                      'krmax', [1, kr3D(1-srw)], ...
                      'residual', [srw, src]);
                    
pressure = 10 * mega * Pascal;
[S, Smax] = deal(0.2, 0.3);

%% reconstruct fine-scale saturations

% capillary fringe fluid
[h_CP, h_max_CP] = upscaledSat2height(S, Smax, Gt, ...
                                      'pcWG', fluidCP.pcWG, ...
                                      'rhoW', fluidCP.rhoW, ...
                                      'rhoG', fluidCP.rhoG, ...
                                      'p', pressure);

[s_CP, smax_CP, seff_CP] = height2finescaleSat(h_CP, h_max_CP, Gt, srw, src, ...
                                               'invPc3D', fluidCP.invPc3D, ...
                                               'rhoW', fluidCP.rhoW(pressure), ...
                                               'rhoG', fluidCP.rhoG(pressure));

% sharp interface fluid
[h_SI, h_max_SI] = upscaledSat2height(S, Smax, Gt, 'resSat', [srw, src]);
[s_SI, smax_SI, seff_SI] = height2finescaleSat(h_SI, h_max_SI, Gt, srw, src);

%% Plot saturation profiles

figure;
clf; subplot(1,2,1);
plot_saturation_profile(s_CP, seff_CP, smax_CP, G);
title('Capillary fringe');

subplot(1,2,2);
plot_saturation_profile(s_SI, seff_SI, smax_SI, G);
title('Sharp interface model');

set(gcf, 'position', [2720 1582 1072 804]);

%% Plot upscaled relperms and capillary pressure

s_input = linspace(0, 1-srw, 100);
[krvalCP, krvalSI, pcapCP, pcapSI] = deal(s_input * 0);
[krvalCPhyst, krvalSIhyst, pcapCPhyst, pcapSIhyst] = deal(s_input * 0);

% compute upscaled relative permeabilities and capillary pressure functions
% for both fluids, and with and without hysteresis
for ix = 1:numel(s_input)
    % without hysteresis
    krvalCP(ix) = fluidCP.krG(s_input(ix), pressure);
    krvalSI(ix) = fluidSI.krG(s_input(ix), pressure);
    pcapCP(ix) = fluidCP.pcWG(s_input(ix), pressure);
    pcapSI(ix) = fluidSI.pcWG(s_input(ix), pressure);
    
    % with hystersis (we use max(0.5, s) as the maximum historical saturation)
    krvalCPhyst(ix) = fluidCP.krG(s_input(ix), pressure, 'sGmax', max(0.5, s_input(ix)));
    krvalSIhyst(ix) = fluidSI.krG(s_input(ix), pressure, 'sGmax', max(0.5, s_input(ix)));
    pcapCPhyst(ix) = fluidCP.pcWG(s_input(ix), pressure, 'sGmax', max(0.5, s_input(ix)));
    pcapSIhyst(ix) = fluidSI.pcWG(s_input(ix), pressure, 'sGmax', max(0.5, s_input(ix)));

end

figure;
clf; subplot(1,2,1);
% plot upscaled permeability curves
plot(s_input, krvalCP, 'r'); hold on
plot(s_input, krvalSI, 'b');
plot(s_input, krvalCPhyst, 'r--');
plot(s_input, krvalSIhyst, 'b--');
title('Upscaled relative permeability curves.');
legend({'sharp interface, no hyst.', 'cap. fringe, no hyst.', ...
        'sharp interface with hyst.', 'cap. fringe with hyst.'}, ...
       'location', 'northwest');
xlabel('CO_2 saturation');
ylabel('relative pemeability');
set(gca, 'fontsize', 14);

% plot upscaled capillary pressure curves
subplot(1,2,2);
plot(s_input, pcapCP/mega, 'r'); hold on
plot(s_input, pcapSI/mega, 'b');
plot(s_input, pcapCPhyst/mega, 'r--');
plot(s_input, pcapSIhyst/mega, 'b--');
legend({'sharp interface, no hyst.', 'cap. fringe, no hyst.', ...
        'sharp interface with hyst.', 'cap. fringe with hyst.'}, ...
       'location', 'northwest');
xlabel('CO_2 saturation');
ylabel('upscaled capillary pressure (MPa)');
set(gca, 'fontsize', 14);

set(gcf, 'position', [2715 1877 1399 509]);

%% Copyright Notice
%
% <html>
% <p><font size="-1">
% Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.
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
