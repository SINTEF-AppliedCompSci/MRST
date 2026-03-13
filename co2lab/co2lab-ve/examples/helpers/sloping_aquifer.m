function [G, states_3D, G_VE, states_VE, sat_VE3D, timing, wellcell, model, ...
          model_VE, schedule, schedule_VE] = sloping_aquifer(varargin)
%Undocumented Utility Function

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

mrstModule add static-modeling % for generating heterogeneity fields
mrstModule add co2lab-mit % for 3D residual saturaiton
mrstModule add ad-blackoil % for GenericBlackOilModel
gravity on

options.cross_sectional = true;
options.heterogeneity = 'none'; % other options are 'gaussian', 'prograding', and 'retrograding'
options.integrate_heterogeneity = true;
options.xres = 100;
options.yres = 30;
options.zres = 10;
options.srw = 0;
options.srg = 0;
options.cap_fringe = false;
options.pc_scaling_fac = 0.4; % only relevant for capillary fringe models
options.skip_3D = false;
options.avg_perm = 400 * milli * darcy;
options.avg_poro = 0.2;
options.slope = 0.04;
options.thickness = 20 * meter; 
options.trap1_size = 20;
options.trap2_size = 10;
options.mass_injection_rate = 20 * kilo * 1e3 / year;
options.inj_period = 1 * year;
options.inj_steps = 12;
options.migr_period = 2 * year;
options.migr_steps = 12;
options.relperm_exponents = [2, 2]; % default: quadratic relperm curves
options.maxTimestepCuts = 14;

options = merge_options(options, varargin{:});


%% Set up a 3D test grid with two structual traps

[xres, yres, zres] = deal(options.xres, options.yres, options.zres);

if options.cross_sectional
    % we only simulate on a vertical slice
    yres = 1; 
end

[xlen, ylen, zlen] = deal(2000 * meter, 600 * meter, options.thickness);

G = make_testgrid([xres, yres, zres], ...      % x, y and z resolution
                  [xlen, ylen, zlen], ...      % x, y and z lengths
                  5, options.slope, 1000, ...  % bend, slope and depth
                  [0.5, 0.3, options.trap1_size;    % trap 1 x-position, width and size
                   0.75, 0.2, options.trap2_size]); % trap 2 x-position, width and size

%% Set up the rest of the simulation input data

% define rock object, one value per gridcell
if ~strcmpi(options.heterogeneity, 'none')

    expo = 1.5;
    scale = diag([xlen, ylen, zlen]);
    scale = scale * diag([1 1 200]);% rock structure
    res = [xres, yres, zres];
    
    if options.cross_sectional
        scale(2,:) = [];
        scale(:,2) = [];
        res(2) = [];
    end

    layersize = prod(res) / res(end);        
        
    if strcmpi(options.heterogeneity, 'gaussian')
        
        structure_length_scale = 2000; 
            
        % compute an anisotropic Gaussian field to represent log-perm
        logperm = GaussianProcessND(res, ...
                                    @(pt) exp(-sqrt(sum((pt * scale).^2, 2))/structure_length_scale));
        
        % compute permeability from its logarithm, and scale it to ensure exact
        % average
        perm = 10.^(logperm/3); % /5
        perm = perm / mean(perm(:)) * options.avg_perm;
        
    elseif strcmpi(options.heterogeneity, 'prograding')
        % define a permeability that is constant per layer, and increases
        % towards the top of the grid
        perm = reshape(ones(layersize, 1) .* fliplr(1:res(end)), [], 1);
        perm = perm .^ expo;
        perm = perm / mean(perm) * options.avg_perm;
        
    elseif strcmpi(options.heterogeneity, 'retrograding')
        % define a permeability that is constant per layer, and decreases
        %towards the top of the grid
        perm = reshape(ones(layersize, 1) .* (1:res(end)), [], 1);
        perm = perm .^ expo;
        perm = perm / mean(perm) * options.avg_perm;
    else
        % unrecgonized heterogeneity type
        error('Unrecognized heterogeneity type');
    end

    % computing computing porosity using Cozeny-Karman relation
    poro = poroFromPerm(perm, options.avg_poro, 1e-10);
    
    rock = struct('perm', perm(:), 'poro', poro(:));
else
    N = G.cells.num;
    rock.poro = repmat(options.avg_poro, N, 1);
    rock.perm = repmat(options.avg_perm, N, 1);
end

% define initial state

rhow = 1050; % density of brine

pfun = @(z) rhow * norm(gravity) * z;

initState.pressure = pfun(G.cells.centroids(:,3)); % initial pressure
initState.s = repmat([1, 0], G.cells.num, 1); % initial saturations

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

% Use function 'initSimpleADIFluid' to make a simple fluid object
fluid = initSimpleADIFluid('phases', 'WG'                     , ...
                           'mu'  , [muw, muco2]               , ...
                           'rho' , [rhow, rhoc]               , ...
                           'pRef', p_ref                      , ...
                           'c'   , [cf_wat, cf_co2]           , ...
                           'cR'  , cf_rock                    , ...
                           'smin', [options.srw, options.srg] , ...
                           'n'   , options.relperm_exponents); 

% initSimpleADIFluid assumes endpoint relperm is 1.  Fix this by redefining
% relperm functions.
krG = @(sg) sg.^options.relperm_exponents(2);  % drainage curve
krGmax = krG(1-options.srw);
fluid.krG = @(sg) min(krGmax, krG(sg));
fluid.krW = @(sw) ( max(sw-options.srw, 0) ./ (1-options.srw)) .^ options.relperm_exponents(1);


if options.srg > 0 || options.srw > 0 % residual saturation present

    srg = max(options.srg, 1e-6); % Ensure not completely zero to avoid corner
    srw = max(options.srw, 1e-6); % case in HystereticRelativePermeability
    
    sgmax = 1 - srw;
    krGimb = @(s) fluid.krG( ...
        sgmax * max((s-srg) ./ (sgmax - srg), 0));

    fluid.krG = {fluid.krG, krGimb};
    fluid.krW = {fluid.krW, fluid.krW};

    fluid.krPts.w = [srw, srw, 1, 1];
    fluid.krPts.g = [0, 0, 1-srw, fluid.krG{1}(1-srw);
                     0, srg, 1-srw, fluid.krG{2}(1-srw)];

    fluid.krHyst = 2; % the presence of the krHyst field will activate 
                      % hysteretic relative permeability in the model
    %fluid.ehystr = {0.1 2 1 0.01};
    fluid.ehystr = {0.1 2 1 0.1};
    
    rock.regions.saturation = ones(G.cells.num, 1);
    rock.regions.imbibition = 2 * ones(G.cells.num, 1);
end

% identify injection cell
well_ix = ceil(xres/8); % towards the left side in x-direction
well_iy = ceil(2*yres/5); % slightly off-center in y direction
well_iz = zres; % bottom in z-direction

wellcell = sub2ind(G.cartDims, well_ix, well_iy, well_iz);

% define injection well
injection_rate = options.mass_injection_rate / fluid.rhoGS;

W = addWell([], G, rock, wellcell, ...
            'refDepth', G.cells.centroids(wellcell, 3), ... % BHP reference depth
            'type', 'rate', ...                             % inject at constant rate
            'val', injection_rate, ...                      % volumetric injection rate
            'comp_i', [0, 1]);                              % inject CO2, not water

% define boundary conditions
bc = open_boundary_conditions(G, pfun, options.cross_sectional);

% define schedule
inj_period = options.inj_period;
inj_steps = options.inj_steps;

if options.cross_sectional
    inj_period = inj_period * 2;
    inj_steps = inj_steps * 2;
end

schedule = simple_injection_migration_schedule(W, bc, inj_period, inj_steps, ...
                                               options.migr_period, options.migr_steps);

% define VE fluid
[G_VE, G] = topSurfaceGrid(G);
rock_VE = averageRock(rock, G_VE);
[C, alpha] = deal(options.pc_scaling_fac, 0.5); % to specify capillary pressure, if we use fringe
relperm_model = 'sharp_interface_simple';
if ~strcmpi(options.heterogeneity, 'none') && options.integrate_heterogeneity
    relperm_model = 'sharp_interface_integrated';
elseif options.cap_fringe
    relperm_model = 'P-scaled table'; 
    %relperm_model = 'S table'; 
end

if iscell(fluid.krW)
    krmax = [fluid.krW{1}(1), fluid.krG{1}(1-options.srw)]; 
    kr3D = fluid.krG{1};
else
    krmax = [fluid.krW(1), fluid.krG(1-options.srw)];  
    kr3D = fluid.krG;
end

fluid_VE = makeVEFluid(G_VE, rock_VE, relperm_model, ...
                       'co2_mu_ref', muco2, ...
                       'wat_mu_ref', muw, ...
                       'co2_rho_ref', rhoc, ...
                       'wat_rho_ref', rhow, ...
                       'co2_rho_pvt', [cf_co2, p_ref], ...
                       'wat_rho_pvt', [cf_wat, p_ref], ...
                       'pvMult_p_ref', p_ref, ...
                       'pvMult_fac', cf_rock, ...
                       'residual', [options.srw, options.srg], ...
                       'krmax', krmax, ...        % only relevant if sharp interface
                       'invPc3D', [C, alpha], ... % only relevant if cap. fringe 
                       'kr3D', kr3D ...           % only relevant if cap. fringe
                       );
                                                

%% Simulate injection in 3D
if options.cap_fringe
    fluid.pcWG = @(sG) fluid_VE.pc3D(1-sG);
end

if ~options.skip_3D
        
    model = TwoPhaseWaterGasModel(G, rock, fluid, 0, 0);
    model.verbose=true;
    nls = NonLinearSolver('maxTimestepCuts', options.maxTimestepCuts, ...
                          'timeStepSelector', IterationCountTimeStepSelector('targetIterationCount', 12));
    %'useRelaxation', true, ...
    tic;
    [wellSol3D, states_3D] = simulateScheduleAD(initState, model, schedule, ...
                                                'NonLinearSolver', nls);
    timing.sim3D = toc;
    states_3D = [{initState}; states_3D];
else
    states_3D = {};
    model = [];
end

%return;

%% Simulate injection in VE
model_VE = CO2VEBlackOilTypeModel(G_VE, rock_VE, fluid_VE);

initState_VE.pressure = pfun(G_VE.cells.z);
initState_VE.s = repmat([1, 0], G_VE.cells.num, 1);
initState_VE.sGmax = initState_VE.s(:,2);

bc_VE = open_boundary_conditions(G_VE, pfun, options.cross_sectional);

W_VE = convertwellsVE(W, G, G_VE, rock_VE);

schedule_VE = simple_injection_migration_schedule(...
    W_VE, bc_VE, inj_period, inj_steps, options.migr_period, options.migr_steps);

model_VE.verbose=true;

tic;
%nls = NonLinearSolver('useRelaxation', true);
[wellSol_VE, states_VE] = simulateScheduleAD(initState_VE, model_VE, schedule_VE);
                                             
timing.simVE = toc;
states_VE = [{initState_VE}; states_VE];

sat_VE3D = {};
invpc = []; % if reconstructing saturation in a sharp interface setting, do
            % not provide an inverse capillary pressure function.
if options.cap_fringe
    invpc = fluid_VE.invPc3D;
end

for i = 1:numel(states_VE)
    s = states_VE{i}.s(:,2);
    smax = states_VE{i}.sGmax;
    p = states_VE{i}.pressure;
    [h, h_max] = upscaledSat2height(s, smax, G_VE, ...
                                    'resSat', [0, 0], ...
                                    'pcWG', fluid_VE.pcWG, ...
                                    'rhoW', fluid_VE.rhoW, ...
                                    'rhoG', fluid_VE.rhoG, 'p', p);
    sat_VE3D = [sat_VE3D, {height2finescaleSat(h, h_max, G_VE, ...
                                               fluid_VE.res_water, fluid_VE.res_gas, ...
                                               'invPc3D', invpc, ...
                                               'rhoW', fluid_VE.rhoW(p), ...
                                               'rhoG', fluid_VE.rhoG(p))}];
end
