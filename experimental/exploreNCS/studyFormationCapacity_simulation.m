% This script is a combination of existing scripts and routines in
% mrst-co2lab, for the purpose of analyzing a user-specified formation, and
% estimating its CO2 storage capacity.

% If using formations found in Barents Sea, this script assumes you have
% downloaded the datafiles and placed them under
% mrst-core/examples/data/CO2Atlas. Format of filename is as follows:
% "Name_thickness" and "Name_top", where Name is the formation name
% (however not followed by 'fm').

moduleCheck('co2lab', 'mex',  'libgeometry', 'opm_gridprocessing');

mrstVerbose on

studyParameterRanges = false;
UseInteractiveTrapping = false;


%% 1. Load and process a formation.
N = 5;
name = 'Utsirafm'; %'Bjarmeland'; %'Sandnesfm'; %'Utsirafm'; %'Sandnesfm';  %'Johansenfm';
[grdecl] = getAtlasGrid(name, 'coarsening',N);
G = mprocessGRDECL(grdecl{1});  % use mprocessGRDECL if its big.  % requires mex module
G = mcomputeGeometry(G);        % function requires libgeometry module

[Gt, rock] = getFormationTopGrid(name, N);
rock.poro = unique(rock.poro); % a single value. Use full porosity map instead?
if isnan(rock.poro)
    fprintf('Porosity required: Update getAtlasGrid() to include porosity of %s formation.\n', name);
    return
end

% Note, one could also use the following function, however this function
% does not return rock data as output:
% Gt = topSurfaceGrid(G);


%% 2. -- 5.
% see 'studyFormationCapacity.m'


%% *************************************************************************
%% 6. Run a VE simulation, and obtain CO2 inventory breakdown:
% a) compare with and without dissolution
% b) compare bdry types (open, closed, semi)

% Simulation options:

% Continue with interactiveTrapping() and simply run a simulation scenario
% given the injection well location that's selected. This simulation is
% executed by migrateInjection(). Inside this function, incompTPFA() and
% implicitTransport() are both called.

% Also, runSleipner.m contains calls to solveIncompFlowVE(),
% mtransportVE(), and explicitTransportVE().

% On the other hand, the GUI exploreSimulation() uses a call to the
% function simulateScheduleAD(). Each input must be constructed.

mrstVerbose true


% Basic routine to perform VE simulation, using simulateScheduleAD().
% _________________________________________________________________________
% 1) set up initial state.
water_density = 1000; % kg per m3
initState.pressure  = Gt.cells.z * norm(gravity) * water_density;   % hydrostatic pressure, in Pa=N/m^2
initState.s         = repmat([1 0], Gt.cells.num, 1);               % sat of water is 1, sat of CO2 is 0
initState.sGmax     = initState.s(:,2);                             % max sat of CO2 is initially 0
initState.rs        = 0 * initState.sGmax;                          % initially 0

figure;
plotCellData(Gt, initState.pressure)
title('Initial Pressure (hydrostatic)'); axis off equal tight
hcb = colorbar;
hcb.Label.String = 'Pascals'; set(hcb, 'fontSize', 18)


% OR: get literature data
[rho, mu, sr, sw] = getValuesSPE134891();

% _________________________________________________________________________
% 2) set up schedule (wells, bc, etc.).

% WELLS:
rock2D = rock; % TODO: use por/perm maps
rhoCref = 760 * kilogram / meter ^3; % an (arbitrary) reference density


% Specify well location in ED50 coord system (i.e., Sleipner: (x,y) =
% (4.38516e5, 6.47121e6) approximately (Singh et al. 2010)):
wellXcoord = 4.38516e5;
wellYcoord = 6.47121e6;

% get cell index of specified coordinate location:
wellIndexI_candidates = find(Gt.cells.centroids(:,1)<wellXcoord);
wellIndexJ_candidates = find(Gt.cells.centroids(:,2)<wellYcoord);

PossibleIndexI = zeros(Gt.cells.num,1);
PossibleIndexI(wellIndexI_candidates) = Gt.cells.centroids(wellIndexI_candidates,1);

PossibleIndexJ = zeros(Gt.cells.num,1);
PossibleIndexJ(wellIndexJ_candidates) = Gt.cells.centroids(wellIndexJ_candidates,2);

% get the matching PossibleIndexI and PossibleIndexJ that are non-zero, and
% get the cell index of the maximum matching value
matchingIndex = PossibleIndexI.*PossibleIndexJ;
maxMatchingIndex = find(matchingIndex,1,'last');

wellCellIndex = maxMatchingIndex;


% OR: compute distances between all cell centroids to specified well
% location, and use the cell centroid which has the minimum squared norm.
% TODO.
%dv = bsxfun(@minus, Gt.cells.centroids(:,1:2), [wellXcoord, wellYcoord]);
%[v,i] = min(sum(dv.^2, 2));



[i, j] = ind2sub(Gt.cartDims, wellCellIndex);

% Check coordinate the wellCellIndex corresponds to:
wellCoord_x = Gt.cells.centroids(wellCellIndex,1);
wellCoord_y = Gt.cells.centroids(wellCellIndex,2);
wellCoord_z = 0;

annual_inj_rate = 10; % Mt/year
inj_rate = annual_inj_rate * 1e9 / rhoCref / (365*24*60*60); % m^3/s

W_on  = addWell([], Gt.parent, rock2D, wellCellIndex, ...
    'name', sprintf('W%i', i), 'Type', 'rate', 'Val', inj_rate, 'comp_i', [0 1]);
W_off = W_on;
W_off.val = 0;

% Put into schedule fields --> [injection period; migration period]
schedule.control(1).W = W_on;
schedule.control(2).W = W_off;


% BOUNDARY CONDITIONS:
% First get the faces of the boundaries. face.neighbors are the indices of
% the cells on either side of the faces, i.e., face.neighbor(100,1) and
% face.neighbor(100,2) give the index of the cells on either side of face
% with index 100. Any 0 cell index means there is no cell, i.e., the face
% is along an external boundary of the domain. Thus bdryFaces may be
% obtained by finding all the face indices that contain a 0 cell index on
% either side.
bdryFaces = find( Gt.faces.neighbors(:,1).*Gt.faces.neighbors(:,2) == 0 );
% Then use function bc = addBC(bc, faces, type, value, varargin)
bc = addBC( [], bdryFaces, ...
    'pressure', Gt.faces.z(bdryFaces) * water_density * norm(gravity), ...
    'sat', [1 0] );
% Put into schedule fields --> [injection period; migration period]
schedule.control(1).bc = bc;
schedule.control(2).bc = bc;
             

% TIME STEP:
% Specify and compute time step size for injection period.
inj_time = 50 * year; % (if unit not specified, schedule.step.control is treated as seconds)
inj_steps = 10;
dTi = inj_time / inj_steps;

% Specify and compute time step size for migration period. 
mig_time = 3000 * year;
mig_steps = 30;
dTm = mig_time / mig_steps;

% For simulation schedule
istepvec = ones(inj_steps, 1) * dTi;
mstepvec = ones(mig_steps, 1) * dTm;

schedule.step.val       = [istepvec; mstepvec];
schedule.step.control   = [ones(inj_steps, 1); ones(mig_steps, 1) * 2];


% _________________________________________________________________________
% 3) set up model (grid, rock and fluid properties).
seafloor_temp = 7; % Celsius
seafloor_depth = 100; % meters
temp_gradient = 35.6; % Celsius / km
caprock_temperature =  273.15 + seafloor_temp + (Gt.cells.z - seafloor_depth) / 1e3 * temp_gradient; % Kelvin

% water density, pvt info:
water_compr_val = 0; %4.3e-5/barsa; % will convert to compr/Pa
% pressure at which water density equals the reference density:
ref_p           = mean(initState.pressure); % use mean pressure as ref for linear compressibilities
% pore volume multiplier:
pvMult          = 0; %1e-5/barsa;
% residuals:
water_residual  = 0.11;
co2_residual    = 0.21;
% dissolution:
dis_max         = (53 * kilogram / meter^3) / rhoCref; % from CO2store

fluid = makeVEFluid(Gt, rock2D, 'sharp interface', ...
                              'fixedT'      , caprock_temperature, ...
                              'wat_rho_pvt' , [water_compr_val, ref_p], ...
                              'wat_rho_ref' , water_density, ...
                              'pvMult_p_ref', ref_p, 'pvMult_fac', pvMult, ...                            , ...
                              'residual'    , [water_residual,  co2_residual] , ...
                              'dissolution' , false); %, 'dis_max', dis_max);
model = CO2VEBlackOilTypeModel(Gt, rock2D, fluid);


% _________________________________________________________________________
% 4) call to simulateScheduleAD().
[wellSols, states, sim_report] = simulateScheduleAD(initState, model, schedule);


% save this simulation scenario:
% wellSols, states, sim_report, Gt, model

% _________________________________________________________________________
% 5) Look at results:

% BHP VS TIME
time = sim_report.ReservoirTime;
bhp = zeros(numel(wellSols),1);
for i = 1:numel(wellSols)
    bhp(i) = wellSols{i}.bhp; % bhp is in Pa=N/m^2
end
figure;
plot(time/365/24/60/60,bhp,'x--')
xlabel('Reservoir time, years'); ylabel('well bhp, Pascals=10^{-5}bars');


% ACCUM CO2 VS TIME
accumCO2sat = zeros(numel(states),1);
accumCO2mass = zeros(numel(states),1);
for i = 1:numel(states)
    accumCO2sat(i) = sum( states{i}.s(:,2).*model.G.cells.volumes ); % sat vals are in terms of pore volume
    
    satCO2          = states{i}.s(:,2);
    densityCO2      = fluid.rhoG(states{i}.pressure); 
    accumCO2mass(i) = sum( model.rock.poro * model.G.cells.volumes .* model.G.cells.H .* satCO2 .* densityCO2 );
end
figure;
plot(time/365/24/60/60,accumCO2mass,'o-')
xlabel('Reservoir time, years'); ylabel('Accumlated CO2 mass, kg');



% PROFILES AT SELECT TIME
times2plot = [sim_report.ReservoirTime(1) sim_report.ReservoirTime(round(end/2)) sim_report.ReservoirTime(end)];

% meaningful profiles
press       = states{end}.pressure;
pressDiffFromHydrostatic = press - initState.pressure;
densityCO2  = fluid.rhoG(states{end}.pressure);  % fluid.rhoG is function handle to get CO2 density
satCO2      = states{end}.s(:,2);
massCO2     = model.rock.poro*model.G.cells.volumes.* model.G.cells.H.*satCO2.*densityCO2; % kg

bf = boundaryFaces(Gt);

h = figure;
subplot(1,6,1)
plotCellData(Gt, caprock_temperature, 'EdgeColor','none')
title('Caprock Temperature'); axis off equal tight
hcb = colorbar;
hcb.Label.String = 'Kelvin'; set(hcb, 'fontSize', 18)

subplot(1,6,2)
plotCellData(Gt, press, 'EdgeColor','none')
title('Pressure'); axis off equal tight
hcb = colorbar;
hcb.Label.String = 'Pascals'; set(hcb, 'fontSize', 18)

subplot(1,6,3)
plotCellData(Gt, pressDiffFromHydrostatic, 'EdgeColor','none')
title({'Pressure diff. from hydrostatic';'i.e., the initial condition'}); axis off equal tight
hcb = colorbar;
hcb.Label.String = 'Pascals'; set(hcb, 'fontSize', 18)

subplot(1,6,4)
plotCellData(Gt, densityCO2, 'EdgeColor','none')
title('CO2 density at Caprock'); axis off equal tight
hcb = colorbar;
hcb.Label.String = 'kg/m^3'; set(hcb, 'fontSize', 18)

subplot(1,6,5)
%plotGrid(Gt, 'FaceColor', 'white', 'EdgeAlpha', 0)
plotFaces(Gt, bf, 'EdgeColor','k', 'LineWidth',3);
plotCellData(Gt, satCO2, satCO2>(0.1/100), 'EdgeColor','none') %satCO2~=0)
title('Saturation of CO2'); axis off equal tight
hcb = colorbar;
hcb.Label.String = 'Pore Volume CO2 Saturation'; set(hcb, 'fontSize', 18)

subplot(1,6,6)
%plotGrid(Gt, 'FaceColor', 'white', 'EdgeAlpha', 0)
plotFaces(Gt, bf, 'EdgeColor','k', 'LineWidth',3);
plotCellData(Gt, massCO2/1e9, satCO2>(0.1/100), 'EdgeColor','none') % only plot plume that has sat > 0.1 percent
title('Mass of CO2'); axis off equal tight
hcb = colorbar;
hcb.Label.String = 'Mt'; set(hcb, 'fontSize', 18)


%% 7. Use objective function to get optimized injection rates at wells.


