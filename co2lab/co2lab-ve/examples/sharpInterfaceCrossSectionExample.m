mrstModule add co2lab-common co2lab-ve
mrstModule add static-modeling
mrstModule add ad-core ad-props ad-blackoil
mrstModule add mrst-gui
gravity on;

% In this example, we compare the result of a full 3D simulation and its
% corresponding VE simulation on a cross section of a sloping aquifer.  
% The simulation can be run with or without the following effects:
% - Residual trapping
% - Heterogeneous rock properties
% - Impact of small-scale oscillations in caprock ('rugosity')
% 
% Note that if you include sub-scale oscillations in caprock, the grid used 
% in the 3D simulation will have to be much finer in order to properly capture
% the geometry.  The VE model will still use the original resolution grid, as
% caprock rugosity is modeled as an upscaled effect anyway.
%
% You may turn on and off these effects in the 'User options' section below.

%% User options
heterogeneity = true;
residual = true;
rugosity = false;

show_plots = true;

%% Set up 3D test grid

[xres, zres] = deal(100, 25); 
[L, H] = deal(5000 * meter, 25 * meter);

% The call to 'make_testgrid' sets up a grid with the specified dimension,
% resolution and slope, and with a large structural trap towards the right side.
G = make_testgrid([xres, 1, zres], [L, 1, H], 5, 0.02, 2000, [0.6, 900, 30]);

[xfac, GFine] = deal(1, []);
if rugosity
    % we need a high-resolution grid to capture the high-frequency caprock oscillations
    xfac = 12;
    GFine = make_testgrid([xres*xfac, 1, zres], [L, 1, H], 5, 0.02, 2000, [0.6, 900, 30], ...
                          'ripples_amplitude_x', H/30, 'ripples_number_x', 200);
end


%% Set up rock properties

avg_poro = 0.2;
avg_perm = 200 * milli * darcy;

if heterogeneity
    scale = diag([L/xres, H/zres]);
    scale = scale * diag([0.1, 5]); % anisotropic rock structure

    % Compute an anisotropic gaussian field to represent log-perm
    logperm = GaussianProcessND([xres, zres], @(xy) exp(-sqrt(sum((xy * scale).^2, 2))/0.3));
    
    % compute permeability from its logarithm, and scale it to ensure exact average
    perm = 10.^(logperm/5);
    perm = perm / mean(perm(:)) * avg_perm; 

    % compute corresponding porosity using Cozeny-Karman relation
    poro = poroFromPerm(perm, avg_poro, 1e-10);
    
else
    % Rock properties will be constant, and equal to the prescribed averages.
    poro = avg_poro * ones(G.cells.num, 1);
    perm = avg_perm * ones(G.cells.num, 1);
end

% We make two versions of the rock object, one standard, and one high-res.  The 
% latter one will have to be used for the 3D simulation when caprock rugosity is
% turned on.
rock = struct('poro', poro(:), 'perm', perm(:));
rockFine = struct('poro', reshape(repmat(poro(:), 1, xfac)', [], 1), ...
                  'perm', reshape(repmat(perm(:), 1, xfac)', [], 1));

% choosing grid and rock to use for 3D simulation.
Gsim = ifelse(rugosity, GFine, G);
rocksim = ifelse(rugosity, rockFine, rock);

if show_plots
    % Show the cross-section grid with the porosity and permeability fields.
    figure; 
    subplot(2, 1, 1); 
    plotCellData(Gsim, rocksim.poro, 'edgealpha', 0.2); view(0,0); colorbar; title('porosity');
    subplot(2, 1, 2); 
    plotCellData(Gsim, rocksim.perm/darcy, 'edgealpha', 0.2); view(0,0); colorbar; 
    title('permeability (D)');
    
    set(gcf, 'position', [900, 800, 1000, 500]);
end

%% Simulate 3D injection and migration

% Define initial state
rhow = 1050; % density of brine
pfun = @(z) rhow * norm(gravity) * z; % hydrostatic pressure

initState.pressure = pfun(Gsim.cells.centroids(:,3)); % initial pressure
initState.s = repmat([1, 0], Gsim.cells.num, 1); % initial saturations

% define fluid object
co2     = CO2props(); % load sampled tables of co2 fluid properties
p_ref   = mean(initState.pressure); % choose reference pressure
t_ref   = 70 + 273.15; % choose reference temperature, in Kelvin
rhoc    = co2.rho(p_ref, t_ref); % co2 density at ref. press/temp
cf_co2  = 0; % co2 compressibility (zero) 
             % co2.rhoDP(p_ref, t_ref) / rhoc; % co2 compressibility
cf_wat  = 0; % brine compressibility (zero)
cf_rock = 4.35e-5 / barsa; % rock compressibility
muw     = 8e-4 * Pascal * second; % brine viscosity
muco2   = co2.mu(p_ref, t_ref) * Pascal * second; % co2 viscosity
smin = ifelse(residual, [0.2, 0.2], [0, 0]);
fluid = initSimpleADIFluid('phases', 'WG'             , ...
                           'mu'    , [muw, muco2]     , ...
                           'rho'   , [rhow, rhoc]     , ...
                           'pRef'  , p_ref            , ...
                           'c'     , [cf_wat, cf_co2] , ...
                           'cR'    , cf_rock          , ...
                           'smin'  , smin             , ...
                           'n'     , [2 2]); % quadratic relperm curves

% define injection well
well_ix = 1; %ceil(xres/8);
well_iz = zres;
wellcell = sub2ind(Gsim.cartDims, well_ix, 1, well_iz);
injection_rate = 3.0 * kilo * 1e3 / year/ fluid.rhoGS;

W = addWell([], Gsim, rocksim, wellcell, ...
            'refDepth', Gsim.cells.centroids(wellcell, 3), ...
            'type', 'rate', ...
            'val', injection_rate, ...
            'comp_i', [0, 1]);

% define open boundary on right side
open_faces = find(Gsim.faces.centroids(:,1) == max(Gsim.faces.centroids(:,1)));
bc_cells = sum(Gsim.faces.neighbors(open_faces, :), 2);

bc = addBC([], open_faces, 'pressure', initState.pressure(bc_cells), 'sat', [1, 0]);

% define injection schedule (1 year of injection, 100 years of migration)
inj_period = 1 * year;
inj_steps = ifelse(rugosity, 36*5, 36);
migr_period = 100 * year;
migr_steps = ifelse(rugosity, 50*5, 50);

schedule.control = struct('W', W, 'bc', bc);
schedule.control(2) = struct('W', W, 'bc', bc);
schedule.control(2).W.val = 0;

dT_injection = rampupTimesteps(inj_period, inj_period/inj_steps, 7);
dT_migration = repmat(migr_period / migr_steps, migr_steps, 1);

schedule.step.val = [dT_injection; dT_migration];
schedule.step.control = [ones(numel(dT_injection), 1); ...
                         2 * ones(numel(dT_migration), 1)];
                    

% Run simulation
model = TwoPhaseWaterGasModel(Gsim, rocksim, fluid, 0, 0, 'verbose', false);
[wellSol3D, states3D] = simulateScheduleAD(initState, model, schedule);

if show_plots
    % We insert the permeability and porosity fields into the simulation results,
    % to allow them to be visualized as well in the call to 'plotToolbar' below.
    for i=1:numel(states3D)
        states3D{i}.perm = rocksim.perm/darcy;
        states3D{i}.poro = rocksim.poro;
    end
    
    % plotToolbar will open an interactive window, where the result of the simulation
    % can be explored in detail.
    figure; plotToolbar(Gsim, states3D, 'field', 's:2'); view(0,0); title('3D simulation');
end


%% Run VE simulation

% We now seek to simulate the same scenario, this time using the
% vertical-equilibrium model (which is much more computationally efficient).  To
% do this, we will have to convert all the objects (grid, rock fluid, initial
% state, wells and boundary conditions) to their corresponding VE form.  

% Convert grid and rock to VE form
[Gt, G] = topSurfaceGrid(G);
rockVE = averageRock(rock, Gt);

% Establish the VE form of the fluid to use.  The upscaled relperm and capillary
% pressure curves will be modeled depending on the choice of simulation model
% type.  The 'rugosity' variable allows and upscaled representation of the effect
% of caprock oscillations, without having to explicitly model them on a fine grid.
model_type = ifelse(heterogeneity, 'sharp_interface_integrated', ...
                                   'sharp_interface_simple');
dh = ifelse(rugosity, H/30, 0);
fluidVE = makeVEFluid(Gt, rockVE, model_type, ...
                      'residual', smin, ...
                      'co2_mu_ref', muco2, ...
                      'wat_mu_ref', muw, ...
                      'co2_rho_ref', rhoc, ...
                      'wat_rho_ref', rhow, ...
                      'co2_rho_pvt', [cf_co2, p_ref], ...
                      'wat_rho_pvt', [cf_wat, p_ref], ...
                      'pvMult_p_ref', p_ref, ...
                      'rugosity', dh, ...
                      'pvMult_fac', cf_rock);
                      %'krmax', (1 - fliplr(smin)).^1); % corey exponent 2

% Here we create the simulation model.
modelVE = CO2VEBlackOilTypeModel(Gt, rockVE, fluidVE);

% The VE form of the initial state ressembles that of the 3D form, but 
% has a lot less cells (since the z-dimension is integrated out).
initStateVE.pressure = pfun(Gt.cells.z);
initStateVE.s = repmat([1, 0], Gt.cells.num, 1);
initStateVE.sGmax = initStateVE.s(:,2);

% Define VE boundary conditions (open boundary on right side), in the same way
% as for the 3D c ase, but using the VE grid.
open_faces_VE = find(Gt.faces.centroids(:,1) == max(Gt.faces.centroids(:,1)));
bc_cells_VE = sum(Gt.faces.neighbors(open_faces_VE,:), 2);
bcVE = addBC([], open_faces_VE, ...
             'pressure', initStateVE.pressure(bc_cells_VE), ...
             'sat', [1, 0]);

% Since the 3D simulation was run on a finer 3D grid if rugosity was included, we
% cannot simply convert the existing well, but have to generate a corresponding
% well on the original 3D grid and then convert it to VE.  If rugosity was turned
% off, we can convert the well directly.
if rugosity
    wellcell_tmp = sub2ind(G.cartDims, well_ix, 1, well_iz);
    Wtmp = addWell([], G, rock, wellcell_tmp, ...
            'refDepth', G.cells.centroids(wellcell_tmp, 3), ...
            'type', 'rate', ...
            'val', injection_rate, ...
            'comp_i', [0, 1]);
    WVE = convertwellsVE(Wtmp, G, Gt, rockVE);
else
    % we can recycle the existing well, as it was defined on the parent grid
    WVE = convertwellsVE(W, G, Gt, rockVE);
end

% We use the same schedule as before, but with the VE versions of wells and
% boundary conditions.
scheduleVE.control = struct('W', WVE, 'bc', bcVE);
scheduleVE.control(2) = struct('W', WVE, 'bc', bcVE);
scheduleVE.control(2).W.val = 0;
scheduleVE.step = schedule.step;
                      
% run VE simulation
[wellSolVE, statesVE] = simulateScheduleAD(initStateVE, modelVE, scheduleVE);

% Reconstruct 3D solution from VE solution (so we can plot it on the original 3D grid).
for i = 1:numel(statesVE)
    s = statesVE{i}.s(:,2);
    smax = statesVE{i}.sGmax;
    [h, hmax] = upscaledSat2height(s, smax, Gt, 'resSat', smin); 
    statesVE{i}.sg3D = height2finescaleSat(h, hmax, Gt, smin(1), smin(2));
    statesVE{i}.perm = rock.perm/darcy;
    statesVE{i}.poro = rock.poro;
end

if show_plots
    % Open a new window to explore the result of the VE simulation.  If you left
    % open the 'plotToolbar' window created above for the 3D case, you can now easily
    % compare the two simulation outcomes.
    figure; plotToolbar(G, statesVE, 'field', 'sg3D'); view(0,0); title('VE simulation');
end

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
