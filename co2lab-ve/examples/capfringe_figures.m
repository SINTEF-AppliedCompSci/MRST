% Plot graphs of the following, with and without residual trapping:
% 
% - saturation vs. depth
% - relperm as fct. of upscaled saturation
% - capillary pressure as fct. of upscaled saturation

%% Make grid

% make a single column grid as a basis for the test; use a very high vertical
% resolution (10 cm thick cells)
G = cartGrid([1, 1, 1000], [1, 1, 100] * meter);
[Gt, G] = topSurfaceGrid(G);

%% Create fluids object with upscaled saturation, relperm and capillary
%% pressure functions

srw = 0.2; %0.2;
src = 0.25;
rock2D.poro = 0.3;
T = 300;
[C, alpha] = deal(0.4, 0.5); % specifies the capillary pressue function
kr3D = @(s) max((s-src) ./ (1-src), 0) .^ 2;

% fluid with capillary fringe
fluidCP = makeVEFluid(Gt, rock2D, 'P-scaled table', ...
                    'fixedT', T, ...
                    'residual', [srw, src], ...
                    'invPc3D', [C, alpha], ...
                    'kr3D', kr3D);

% simple, sharp-interface fluid
fluidSI = makeVEFluid(Gt, rock2D, 'sharp_interface_simple', ...
                      'fixedT', T, ...
                      'krmax', [1, kr3D(1-srw)], ...
                      'residual', [srw, src]);
                    
pressure = 10 * mega * Pascal;
[S, Smax] = deal(0.1, 0.15);

%% reconstruct fine-scale saturations

% capillary fringe fluid
[h_CP, h_max_CP] = upscaledSat2height(S, Smax, Gt, ...
                                      'pcWG', fluidCP.pcWG, ...
                                      'rhoW', fluidCP.rhoW, ...
                                      'rhoG', fluidCP.rhoG, ...
                                      'p', pressure);

[s_CP, smax_CP, seff_CP] = height2Sat(h_CP, h_max_CP, Gt, srw, src, ...
                                      'invPc3D', fluidCP.invPc3D, ...
                                      'rhoW', fluidCP.rhoW(pressure), ...
                                      'rhoG', fluidCP.rhoG(pressure));

% sharp interface fluid
[h_SI, h_max_SI] = upscaledSat2height(S, Smax, Gt, 'resSat', [srw, src]);
[s_SI, smax_SI, seff_SI] = height2Sat(h_SI, h_max_SI, Gt, srw, src);

%% Plot saturation profiles

figure;
clf; subplot(1,2,1);
plot_saturation_profile(s_CP, seff_CP, smax_CP, G);

subplot(1,2,2);
plot_saturation_profile(s_SI, seff_SI, smax_SI, G);

%% Plot upscaled relperms and capillary pressure

s_input = linspace(0, 1-srw, 100);
[krvalCP, krvalSI, pcapCP, pcapSI] = deal(s_input * 0);
[krvalCPhyst, krvalSIhyst, pcapCPhyst, pcapSIhyst] = deal(s_input * 0);

for ix = 1:numel(s_input)
    % without hysteresis
    krvalCP(ix) = fluidCP.krG(s_input(ix), pressure);
    krvalSI(ix) = fluidSI.krG(s_input(ix), pressure);
    pcapCP(ix) = fluidCP.pcWG(s_input(ix), pressure);
    pcapSI(ix) = fluidSI.pcWG(s_input(ix), pressure);
    
    % with hystersis
    krvalCPhyst(ix) = fluidCP.krG(s_input(ix), pressure, 'sGmax', max(0.5, s_input(ix)));
    krvalSIhyst(ix) = fluidSI.krG(s_input(ix), pressure, 'sGmax', max(0.5, s_input(ix)));
    pcapCPhyst(ix) = fluidCP.pcWG(s_input(ix), pressure, 'sGmax', max(0.5, s_input(ix)));
    pcapSIhyst(ix) = fluidSI.pcWG(s_input(ix), pressure, 'sGmax', max(0.5, s_input(ix)));

end

figure;
plot(s_input, krvalCP, 'r'); hold on
plot(s_input, krvalSI, 'b');
plot(s_input, krvalCPhyst, 'r--');
plot(s_input, krvalSIhyst, 'b--');

figure;
plot(s_input, pcapCP, 'r'); hold on
plot(s_input, pcapSI, 'b');
plot(s_input, pcapCPhyst, 'r--');
plot(s_input, pcapSIhyst, 'b--');

