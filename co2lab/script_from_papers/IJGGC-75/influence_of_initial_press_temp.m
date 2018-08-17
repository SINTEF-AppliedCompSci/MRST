%% Influence of pressure and temperature on structural trapping estimates
% The following script investigates the influence of initial pressure and
% temperature on the static estimation (i.e., without simulation) of the
% structural trapping capacity in the southern part of Utsira. First,
% structural capacity is calculated for the base model, which assumes
% hydrostatic pressure and a caprock temperature calculated using a thermal
% gradient of 35.6 C/km. Then, pressure is deviated between -20% to 20% of
% its pore-volume weighted average. And then, temperature is deviated by
% considering thermal gradients between 34 and 43 C/km. Plots are generated
% to show how structural trapping capacity changes with initial pressure
% and temperature. Map plots are also generated to show structural trapping
% capacity in the following cases:
%   * base case,
%   * given a pressure deviation of +15%,
%   * given a thermal gradient of 34 C/km.

mrstModule add co2lab
moduleCheck('ad-core','ad-props')


%% Base model case
% assuming:
%   * hydrostatic pressure
%   * temperature assuming temp. gradient of 35.6 C/km

[Gt_base, rock_base, seainfo] = makeMyGeomodel('modify_base_rock',false);
ta_base = trapAnalysis(Gt_base, true); % independent of rock

surface_pressure = 1*atm;


% Plot initial conditions of base model:

% initial pressure
p_init = (seainfo.water_density * norm(gravity()) * Gt_base.cells.z + surface_pressure);
figure, title('Initial Pressure (bars)')
plotCellData(Gt_base, convertTo(p_init,barsa), 'edgecolor','none');
colorbar; axis equal tight off

% initial temperature
t_init = seainfo.seafloor_temp + (Gt_base.cells.z - seainfo.seafloor_depth) ...
    .* seainfo.temp_gradient ./ 1000; % Celsius
t_init = t_init + 273.14; % Kelvin
figure, title('Initial Temperature (C)')
plotCellData(Gt_base, t_init-273.14, 'edgecolor','none');
colorbar; axis equal tight off


% Map of corresponding structural traps and their capacities. Red shaded
% region indicates where CO2 is non-dense according according to (p,t)
% conditions.
p = p_init;
t = t_init;
[~, strap, ~, ~, ~] = compute_trapcap(Gt_base, ta_base, rock_base, ...
    seainfo, surface_pressure);

plotTrapCapacities(Gt_base, strap, p, t);


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
tot_strap = zeros(1,numel(press_deviation_all));
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
% +15% of its pore-volume weighted average (i.e., +11.5 bars). Red shaded
% region indicates where CO2 is non-dense according according to (p,t)
% conditions.
p_dev = 15;
p = p_init + ref_p * p_dev/100;
t = t_init;
[~, strap, ~, ~, ~] = compute_trapcap(Gt_base, ta_base, rock_base, ...
    seainfo, surface_pressure, 'press_deviation',p_dev);

plotTrapCapacities(Gt_base, strap, p, t);


%% Influence of initial temperature on structural trapping estimate
% Then we deviate caprock temperature by considering thermal gradients
% ranging from 34 to 43 degrees C/km.

temp_deviation_all = [34:0.2:43] - seainfo.temp_gradient; % C/km


% Loop through realizations (deviate a parameter)
press_deviation = 0;
seainfo_r = seainfo;
tot_strap = zeros(1,numel(temp_deviation_all));
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
% of 34 C/km. Red shaded region indicates where CO2 is non-dense according
% according to (p,t) conditions
temp_grad = 34; % C/km
t = seainfo.seafloor_temp + (Gt_base.cells.z - seainfo.seafloor_depth) ...
    .* temp_grad ./ 1000; % Celsius
t = t + 273.14; % Kelvin
p = p_init;
[~, strap, ~, ~, ~] = compute_trapcap(Gt_base, ta_base, rock_base, ...
    seainfo, surface_pressure, 'caprockTemp',t);

plotTrapCapacities(Gt_base, strap, p, t);