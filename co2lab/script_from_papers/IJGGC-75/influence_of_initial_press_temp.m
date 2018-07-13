%% Influence of pressure and temperature on structural trapping estimates
% The following script investigates the influence of initial pressure and
% temperature on the static estimation (i.e., without simulation) of the
% structural trapping capacity. First, structural capacity is calculated
% for the base model, which assumes hydrostatic pressure and a caprock
% temperature calculated using a thermal gradient of 35.6 C/km. Then,
% pressure is deviated between -20% to 20% of its pore-volume weighted
% average. And then, temperature is deviated by considering thermal
% gradients between 34 and 43 C/km. Plots are generated to show how
% structural trapping capacity changes with initial pressure and
% temperature. Map plots are also generated to show structural trapping
% capacity in the following cases:
%   * base case,
%   * given a pressure deviation of +15%,
%   * given a thermal gradient of 34 C/km.


moduleCheck('ad-core','co2lab','mrst-gui','ad-props')


%% Base model case
% assuming:
%   * hydrostatic pressure
%   * temperature assuming temp_grad of 35.6 C/km

[Gt_base, rock_base, seainfo] = makeMyGeomodel('modify_base_rock',false);
ta_base = trapAnalysis(Gt_base, true); % independent of rock

surface_pressure = 1*atm;


% Plot initial conditions of base model

% initial pressure
p_init = (seainfo.water_density * norm(gravity()) * Gt_base.cells.z + surface_pressure);
figure, title('Initial Pressure (bars)')
plotCellData(Gt_base, convertTo(p_init,barsa), 'edgecolor','none');
colorbar; axis equal tight off

% initial temperature
t = seainfo.seafloor_temp + (Gt_base.cells.z - seainfo.seafloor_depth) ...
    .* seainfo.temp_gradient ./ 1000; % Celsius
t = t + 273.14; % Kelvin
figure, title('Temperature (C)')
plotCellData(Gt_base, t-273.14, 'edgecolor','none');
colorbar; axis equal tight off

% corresponding structural traps and their capacities
[~, strap, ~, ~, ~] = compute_trapcap(Gt_base, ta_base, rock_base, ...
    seainfo, surface_pressure);
figure, title([num2str(sum(strap)/giga), ' Mt'])
plotFaces(Gt_base, boundaryFaces(Gt_base));
plotCellData(Gt_base, strap/giga, strap > 0, 'edgecolor','none');
colorbar; axis equal tight off;
set(gca,'FontSize',14);

% include red shading to indicate where CO2 is non-dense as well as
% non-supercritical
hold on;
p = p_init;
suitable_cells = aquiferConditionSuitability(p, t, 'plot',false);
non_suitable_cells = ~suitable_cells; 
plotCellData(Gt_base, Gt_base.cells.z, non_suitable_cells, ...
    'facecolor','red', 'facealpha',0.2, 'edgecolor','none');
set(gca,'CLim',clim);
H = removeCells(Gt_base.parent, find(~non_suitable_cells));
H = topSurfaceGrid(H);
plotFaces(H, boundaryFaces(H), 'LineWidth',2, 'FaceColor','r', 'EdgeColor','r');



%% Influence of initial pressure on structural trapping estimate

% We deviate hydrostatic pressure uniformly from -20% to 20% of its
% pore-volume weighted average. This corresponds to a deviation of roughly
% -15 to 15 bars since the reference pressure (i.e., pore-volume weighted
% average pressure) for this geomodel is ~80 bars.
if ~isfield(rock_base,'ntg'), ntg = 1; end
pv = rock_base.poro .* Gt_base.cells.volumes .* Gt_base.cells.H .* ntg; % m^3
ref_p = sum(p_init.*pv) / sum(pv);
press_deviation_all = [-20:1:20];


% Loop through realizations (deviate a parameter)
for r = 1:numel(press_deviation_all)

    press_deviation = press_deviation_all(r);
    
    % Get structural trapping capacity only:
    [~, strap, ~, ~, ~] = compute_trapcap(Gt_base, ta_base, rock_base, ...
        seainfo, surface_pressure, 'press_deviation',press_deviation);

    % total:
    tot_strap(r) = sum(strap)/giga; % Mt

end


% Plot results

% Capacity vs. pressure deviation
% Note: pressure deviation plotted here is in terms of a pressure
% difference from the reference pressure (i.e., pore-volume weighted
% average)
figure, plot(tot_strap, convertTo(ref_p,barsa).*press_deviation_all./100, 'x')
xlabel('CO2 (Mt)')
ylabel('Pressure deviation (bars)')


% Map of structural traps and their capacities given pressure deviation of
% +15% of its pore-volume weighted average (i.e., +11.5 bars)
p_dev = 15;
[~, strap, ~, ~, ~] = compute_trapcap(Gt_base, ta_base, rock_base, ...
    seainfo, surface_pressure, 'press_deviation',p_dev);
figure, title([num2str(sum(strap)/giga), ' Mt'])
plotFaces(Gt_base, boundaryFaces(Gt_base));
plotCellData(Gt_base, strap/giga, strap > 0, 'edgecolor','none');
colorbar; axis equal tight off;
set(gca,'FontSize',14); clim = get(gca,'CLim');

% Include red shading to indicate where CO2 is non-dense as well as
% non-supercritical
hold on;
p = p_init + ref_p * p_dev/100;
suitable_cells = aquiferConditionSuitability(p, t, 'plot',false);
non_suitable_cells = ~suitable_cells; 
plotCellData(Gt_base, Gt_base.cells.z, non_suitable_cells, ...
    'facecolor','red', 'facealpha',0.2, 'edgecolor','none');
set(gca,'CLim',clim);
H = removeCells(Gt_base.parent, find(~non_suitable_cells));
H = topSurfaceGrid(H);
plotFaces(H, boundaryFaces(H), 'LineWidth',2, 'FaceColor','r', 'EdgeColor','r');


%% Influence of initial temperature on structural trapping estimate

% Then we deviate caprock temperature by considering thermal gradients
% ranging from 34 to 43 degrees C/km.
temp_deviation_all = [34:0.2:43] - seainfo.temp_gradient; % C/km


% Loop through realizations (deviate a parameter)
press_deviation = 0;
seainfo_r = seainfo;
clear tot_strap
for r = 1:numel(temp_deviation_all)

    seainfo_r.temp_gradient = seainfo.temp_gradient + temp_deviation_all(r);

    % Get structural trapping capacity only:
    [~, strap, ~, ~, ~] = compute_trapcap(Gt_base, ta_base, rock_base, ...
        seainfo_r, surface_pressure, 'press_deviation',press_deviation);

    % total:
    tot_strap(r) = sum(strap)/giga; % Mt

end


% Plot results
figure, plot(tot_strap, temp_deviation_all + seainfo.temp_gradient, 'x')
xlabel('CO2 (Mt)')
ylabel('Temperature gradient (C/km)')


% Map of structural traps and their capacities given temperature gradient
% of 34 C/km
temp_grad = 34; % C/km
t = seainfo.seafloor_temp + (Gt_base.cells.z - seainfo.seafloor_depth) ...
    .* temp_grad ./ 1000; % Celsius
t = t + 273.14; % Kelvin
[~, strap, ~, ~, ~] = compute_trapcap(Gt_base, ta_base, rock_base, ...
    seainfo, surface_pressure, 'caprockTemp',t);
figure, title([num2str(sum(strap)/giga), ' Mt'])
plotFaces(Gt_base, boundaryFaces(Gt_base));
plotCellData(Gt_base, strap/giga, strap > 0, 'edgecolor','none');
colorbar; axis equal tight off;
set(gca,'FontSize',14); clim = get(gca,'CLim');

% Include red shading to indicate where CO2 is non-dense as well as
% non-supercritical
hold on;
suitable_cells = aquiferConditionSuitability(p_init,t, 'plot',false);
non_suitable_cells = ~suitable_cells; 
plotCellData(Gt_base, Gt_base.cells.z, non_suitable_cells, ...
    'facecolor','red', 'facealpha',0.2, 'edgecolor','none');
set(gca,'CLim',clim);
H = removeCells(Gt_base.parent, find(~non_suitable_cells));
H = topSurfaceGrid(H);
plotFaces(H, boundaryFaces(H), 'LineWidth',2, 'FaceColor','r', 'EdgeColor','r');
