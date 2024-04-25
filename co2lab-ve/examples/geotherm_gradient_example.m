%% Geothermal gradient example
% This example shows the influence the geothermal gradient might have on the
% long-term upslope migration of CO2.  We use the Gassum formation as an
% example, a large sloping open aquifer whose depth varies from more than 3000 m
% to less than 500 m.  Accordingly, the aquifer presents a wide range of
% temperature and pressure conditions.  The impact of the thermal variation
% in the aquifer on fluid properties (and hence migration) may not be
% adequately captured by assuming a fixed mean temperature in the aquifer
% during simulation.  In this script, two simulations are run; one that
% assumes a spatially varying temperature field as dictated by the thermal
% gradient, another where a fixed "mean" temperature is assumed for the whole
% aquifer.  
mrstModule add ad-core ad-props co2lab-common co2lab-ve co2lab-spillpoint mrst-gui

gravity reset on;

%% Set up grid and initial pressure and temperature fields

% We load the formation top grid.  We also specify that the corresponding 3D
% (parent) grid should have five vertical layers - we will use this for
% visualizing the 3D reconstruction of the results at the end of the simulation.
[Gt, rock2D] = getFormationTopGrid('Gassumfm', 2, 'vertical_layers', 5);

% Compute reservoir temperature field.  We assume a top surface temperature
% of 10 degrees C, and a thermal gradient of 25 degrees per kilometer.  This
% allows us to compute the temperature of the aquifer top surface, which we
% assume constant over the vertical thickness of the aquifer.
T_surface = 10 + 273.15;      % top surface temperature, in Kelvin
T_grad = 25 / (kilo * meter); % hypothetical thermal gradient
T_reservoir = T_surface + Gt.cells.z * T_grad;

% For comparison, we also define a reservoir thermal field that is constant, and
% equal to the mean temperature of the reservoir.
T_reservoir_const = repmat(mean(T_reservoir), Gt.cells.num, 1);

% Compute hydrostatic reservoir pressure field (equal for both simulations)
P_surface = 1 * atm;
rhoW_approx = 1025 * kilo * gram / meter^3; % approximate brine density (for
                                            % reservoir pressure estimation)
P_reservoir = P_surface + rhoW_approx * norm(gravity) * Gt.cells.z; 

% Plot reservoir temperature field
figure(1); clf; set(gcf, 'position', [100, 100, 1200, 700]);
plotCellData(Gt, T_reservoir - 273.15, 'edgealpha', 0.2); 
colorbar, axis tight
view(-73, 82)
title('Estimated reservoir temperature (C)');

% Plot reservoir pressure field
figure(2); clf; set(gcf, 'position', [300, 300, 1200, 700]);
plotCellData(Gt, P_reservoir / 1e6, 'edgealpha', 0.2); % plot reservoir pressure in MPa
colorbar, axis tight
view(-73, 82)
title('Estimated initial reservoir pressure (MPa)');

%% Setup fluid object

% Specify the range of coverage of sampled property tables.  If these are not 
% already present, they will be constructed using `CoolProp` (this may take
% some time, and you may be prompted to accept installation of CoolProp).

p_range = [0.1, 60] * mega * Pascal; % pressure range for sampled tables
t_range = [  4, 100] + 274;           % temperature range for sampled tables
   
% Specify residual saturations
srw = 0.2; % residual water saturation
src = 0.2; % residual CO2 saturation

% Make two fluid objects: first fluid uses the full temperature field, second 
% uses the constant temperature field
fluids = {}; 
for tfield = {T_reservoir, T_reservoir_const}

    fluid = makeVEFluid(Gt, rock2D, 'sharp_interface_simple', ...
                        'co2_rho_pvt', [p_range, t_range], ...
                        'co2_mu_pvt', [p_range, t_range], ...
                        'wat_rho_pvt', [p_range, t_range], ...
                        'wat_mu_pvt', [p_range, t_range], ...
                        'residual', [srw, src], ...
                        'reservoirT', tfield{:});
    fluids = [fluids, {fluid}];
end

%% Plot CO2 density at initial reservoir conditions

% We compare CO2 density at initial reservoir conditions under the two
% different thermal fields.
co2rho_init = cellfun(@(x) x.rhoG(P_reservoir), fluids, 'uniformoutput', false);

figure(3); clf; set(gcf, 'position', [350, 350, 2100, 700]);
subplot(1,2,1)
plotCellData(Gt, co2rho_init{1}, 'edgealpha', 0.2); colorbar, axis tight
view(-73, 82)
title('CO2 density at reservoir conditions kg/m3 (with thermal gradient)');

subplot(1,2,2)
plotCellData(Gt, co2rho_init{2}, 'edgealpha', 0.2); colorbar, axis tight
view(-73, 82)
title('CO2 density at reservoir conditions kg/m3 (fixed global temperature)');

% When temperature depends on depth, we see a sharp transition between a
% "dense" and a "gas" phase region of CO2, occuring at around 600 m.  The
% transition between gas and supercritical is more gradual in the case with
% constant temperature.  The contrast between the "gas" and "dense" regions
% in the first plot is so strong that it is hard to see the density
% variations within the dense region.  To make this variation more apparent,
% we plot the field again, excluding any values seen around or above the
% transition line:
co2rho_init_modif = co2rho_init{1};
co2rho_init_modif(co2rho_init{1} < 700) = NaN;

figure(4); clf; set(gcf, 'position', [350, 350, 1000, 700]);
plotCellData(Gt, co2rho_init_modif, 'edgealpha', 0.2); colorbar, axis tight
view(-73, 82)
title(['CO2 density at reservoir conditions kg/m3 (with thermal gradient), ' ...
       'dense region only']);

%% Set up simulation schedule

% Identifying well cell to be the closest to the given location.  This
% location is found in the middle of the domain, around the deeper parts of
% the aquifer, and close to a saddle point of the top surface.
wloc = 1e6 * [0.76, 6.389];
wdist = Gt.cells.centroids - wloc;
[~, wcell] = min(sum(wdist.^2, 2)); % pick the closest cell

% Definining the injection rate
co2_rho_ref = fluids{1}.rhoGS; %@@co2rho_init{1}(wcell);
inj_rate = 5 * mega * tonne / year / co2_rho_ref; % A high injection of 5
                                                  % MT/year, converted to volume

% Defining the well.  Since we already use the top surface grid and the 2D rock
% object for this, we do not subsequently have to convert it using
% `convertwellsVE`.
W = addWell([], Gt, rock2D, wcell, ...
            'name', 'injector', ...
            'type', 'rate', ...
            'val', inj_rate, ...
            'comp_i', [0, 1]); 

% boundary conditions (constant hydrostatic pressure)
boundary_faces = find(prod(Gt.faces.neighbors, 2) == 0);
boundary_cells = sum(Gt.faces.neighbors(boundary_faces,:), 2);
bc = addBC([], boundary_faces, 'pressure', P_reservoir(boundary_cells));
bc.sat = repmat([1, 0], numel(boundary_faces), 1);

clear schedule
schedule.control(1) = struct('W', W, 'bc', bc); % injection phase
schedule.control(2) = struct('W', W, 'bc', bc); % migration phase
schedule.control(2).W.val = 0; % no injection during migration phase

schedule.step.val = [repmat(year, 50, 1); ...
                     repmat(10*year, 30, 1)];
schedule.step.control = [ones(50, 1); ones(30, 1) * 2];


%% Setup initial state, simulation models and run both simulations

initState = struct('pressure', P_reservoir, ...
                   's', repmat([1, 0], Gt.cells.num, 1), ...
                   'sGmax', zeros(Gt.cells.num, 1), ...
                   'rs', zeros(Gt.cells.num, 1));

model_Tfield = CO2VEBlackOilTypeModel(Gt, rock2D, fluids{1});
model_Tfixed = CO2VEBlackOilTypeModel(Gt, rock2D, fluids{2});

[wellSol, states_Tfield] = simulateScheduleAD(initState, model_Tfield, schedule);
[wellSol, states_Tfixed] = simulateScheduleAD(initState, model_Tfixed, schedule);

states_Tfield = [{initState}; states_Tfield];
states_Tfixed = [{initState}; states_Tfixed];


%% Inspect result

% Plot end states

% We can here note that the final CO2 distribution is different between the two
% simulations, with CO2 having advanced further into a second trap in the second
% case (fixed global temperature).
figure(5); clf;
subplot(1,2,1); plotCellData(Gt, states_Tfield{end}.s(:,2), 'edgealpha', 0.2); 
title('End state, with thermal gradient');
subplot(1,2, 2); plotCellData(Gt, states_Tfixed{end}.s(:,2), 'edgealpha', 0.2);
title('End state, with fixed global temperature');
set(gcf, 'position', [150, 150, 1600, 700]);

% Compare inventories
ta = trapAnalysis(Gt, false);
reports_Tfield = postprocessStates(Gt, states_Tfield, model_Tfield.rock, fluids{1}, schedule, ta, []);
reports_Tfixed = postprocessStates(Gt, states_Tfixed, model_Tfixed.rock, fluids{2}, schedule, ta, []);
h1 = figure(6); set(h1, 'position', [100 100 1600 700]);
subplot(1,2,1); plot(1); ax = get(h1, 'currentaxes');
plotTrappingDistribution(ax, reports_Tfield, 'legend_location', 'southeast');
title('Trapping inventory (with thermal gradient)');
subplot(1,2,2); plot(1); ax = get(h1, 'currentaxes');
plotTrappingDistribution(ax, reports_Tfixed, 'legend_location', 'southeast');
title('Trapping inventory (fixed global temperature)');


% Interactive inspection: plotToolbar with 3D grid (complete mapping of states from VE to 3D)
states_Tfield_3D = VEstates23D(states_Tfield, Gt, fluids{1}); % convert to 3D
figure(7)
plotToolbar(Gt.parent, states_Tfield_3D); set(gcf, 'position', [100, 100, 1600, 800]);
title('Simulation results, with thermal gradient, 3D reconstruction');

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
