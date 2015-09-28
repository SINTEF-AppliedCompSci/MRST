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
N = 2;
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


%% 2. Visualize formation grid and top surface.
figure;
plotGrid(G, 'EdgeAlpha', 0.1, 'FaceColor', 'none')
plotCellData(Gt, Gt.cells.z) % depth is in units of meters

view(3);
set(gca,'DataAspect',[1 1 1/100]);
hcb = colorbar;
hcb.TickLabels = hcb.Ticks/1000;
hcb.Label.String = 'Depth below sealevel, km';
xlabel('x-coordinate'); ylabel('y-coordinate'); zlabel('z, m'); grid



%% 3. Analyze traps with spill-point analysis, and plot structural traps.
ta = trapAnalysis(Gt, false); % true for cell-based method

plotCellData(Gt, ta.traps, ta.traps~=0, 'FaceColor', 'r', 'EdgeColor','none')
bf = boundaryFaces(Gt, ta.traps~=0);
plotFaces(Gt, bf, 'EdgeColor','k', 'LineWidth',3);

ta_volumes = volumesOfTraps(Gt, ta);

fprintf('Coarsening level %d:\n', N);
fprintf('  Num. global traps: %d\n', numel(ta_volumes));
fprintf('  Total trap volume: %e m3\n', sum(ta_volumes));
fprintf('  Avg. global trap size: %e m3\n', mean(ta_volumes));

% plot of structural traps in red:
title({[grdecl{1}.name ', coarsening level=' num2str(N)]; ...
    [num2str(numel(ta_volumes)) ' structural traps shown in red']}, 'FontSize', 14)

% ****************
% We can use plotCellData and color the faces according to the trap volume
% in km^3. The cell of Gt is associated with a particular trap by ta.traps,
% and is 0 otherwise. The volume of a particular trap is given by
% ta_volumes.
figure;
plotGrid(G, 'EdgeAlpha', 0.1, 'FaceColor', 'none')

trapcells = ta.traps~=0;
cellsTrapVol = zeros(Gt.cells.num,1);
cellsTrapVol(trapcells) = ta_volumes(ta.traps(trapcells));
plotCellData(Gt, cellsTrapVol/1e3/1e3/1e3, cellsTrapVol~=0)

set(gca,'DataAspect',[1 1 1/100])
hcb = colorbar;
hcb.TickLabels = hcb.Ticks;
hcb.Label.String = 'Trap Volume, km^3';
xlabel('x-coordinate'); ylabel('y-coordinate'); zlabel('z, m'); grid

% Note that it may be more useful to color the faces according to the
% trapped mass of CO2. This requires using the CO2 density data which can
% vary with depth (since pressure and temperature gradients exist in the
% formation), to determine the cooresponding CO2 mass trapped by the traps.
% Porosity might also be required, in order to compute both trapping
% volumes and CO2 mass trapping.
% See functions which produce Utsira plots of traps, catchment area, etc:
%   co2 = CO2props();
%   utsiraCapacity('rhoCfun', @co2.rho);
% The above functions also compute breakdown of CO2 trapping types
% (structural, residual, and dissolution) without VE simulation.






%% 3. b) get trapping breakdown (structural, residual, dissoluion)
% The following section computes the breakdown of structural, residual, and
% dissolution trapping in the formation, assuming CO2 is injected and fills
% the formation completely, i.e., all spill paths are utilized. The
% contents of this section closely follows the implementation found in
% exploreCapacity(), which is a GUI. Note: some fluid parameters are fixed
% to those used in exploreCapacity.

% first, use default values to compute theoretical capacity (upper bound):
out_default = exploreParameterRanges(Gt, rock, ta);

% ****************
% We can visualize the structural traps in terms of CO2 mass, not volume.

% Distributed CO2 mass under structural traps: 
cellsTrapCO2Mass = zeros(Gt.cells.num,1);
cellsTrapCO2Mass(trapcells) = out_default.strap_mass_co2(trapcells);

% Cumulative CO2 mass under structural traps:
trapcaps = accumarray(ta.traps(trapcells), out_default.strap_mass_co2(trapcells));
trapcap_tot = zeros(Gt.cells.num,1); %ones(size(ta.traps)) * NaN;
trapcap_tot(trapcells) = trapcaps(ta.traps(trapcells));

% Boundary of formation
bf = boundaryFaces(Gt);

figure;
subplot(1,4,1)
%plotGrid(G, 'EdgeAlpha', 0.1, 'FaceColor', 'none')
plotFaces(Gt, bf, 'EdgeColor','k', 'LineWidth',3);
plotCellData(Gt, cellsTrapCO2Mass/1e9, cellsTrapCO2Mass~=0)

set(gca,'DataAspect',[1 1 1/100])
hcb = colorbar; %hcb.TickLabels = hcb.Ticks;
hcb.Label.String = 'Distributed CO2 Mass under Trap, Mt';
grid; axis tight; set(hcb, 'fontSize', 18); set(gca, 'fontSize', 10);

subplot(1,4,2)
%plotGrid(G, 'EdgeAlpha', 0.1, 'FaceColor', 'none')
plotFaces(Gt, bf, 'EdgeColor','k', 'LineWidth',3);
plotCellData(Gt, trapcap_tot/1e9, trapcap_tot~=0)

set(gca,'DataAspect',[1 1 1/100])
hcb = colorbar; %hcb.TickLabels = hcb.Ticks;
hcb.Label.String = 'Accumulated CO2 Mass under Trap, Mt';
grid; axis tight; set(hcb, 'fontSize', 18); set(gca, 'fontSize', 10);



% ****************
% We can also visualize the reachable structural capacity, which is based
% on the output of trapAnalysis: ta. Then we compute cumulative structural
% trap volume reached when injecting at a given cell
trees = maximizeTrapping(Gt, 'res', ta, 'calculateAll', true, 'removeOverlap', false);
tvols = [trees.value]; %#ok
int_tr = find(ta.trap_regions); %#ok ixs of cells spilling into interior trap
[dummy, reindex] = sort([trees.root], 'ascend'); %#ok

structural_mass_reached = zeros(Gt.cells.num, 1);
for i = 1:numel(ta.trap_z) % loop over each trap

    % ix of cells spilling directly into this trap
    cix = find(ta.trap_regions == i);

    % cell indices of all cells of this trap, and its upstream traps
    aix = find(ismember(ta.traps, [trees(reindex(i)).traps]));

    % compute total structural trap capacity (in mass terms) of these
    % cells, and store result
    structural_mass_reached(cix) = sum(out_default.strap_mass_co2(aix)); %#ok

end

subplot(1,4,3)
plotFaces(Gt, bf, 'EdgeColor','k', 'LineWidth',3);
plotCellData(Gt, structural_mass_reached/1e3/1e6, 'EdgeColor','none');

set(gca,'DataAspect',[1 1 1/100])
hcb = colorbar; %hcb.TickLabels = hcb.Ticks;
hcb.Label.String = 'Reachable structural capacity, Mt';
grid; axis tight; set(hcb, 'fontSize', 18); set(gca, 'fontSize', 10);


% ****************
% We visualize the spill paths between structural traps
subplot(1,4,4)
mapPlot(gcf, Gt, 'traps', ta.traps, 'rivers', ta.cell_lines);
grid; axis equal tight;





% ****************
% Then, use ranges of the input parameters (such as temperature gradient,
% etc...) to compute the range of trapping capacity (such as total and
% structural).

if studyParameterRanges

    % Range of values:
    tgrad_range             = 10:1:50; % degrees C (not Kelvin)
    press_deviation_range   = -50:2:50;
    res_co2_sat_range       = 0:0.01:0.5;
    res_brine_sat_range     = 0:0.01:0.5;
    diss_max_range          = [0:2:100];

    % computes breakdown given range of one parameter (other parameters will be
    % set to their default values)
    [out_tgrad]     = exploreParameterRanges(Gt, rock, ta, 'tgrad',           tgrad_range);
    [out_press]     = exploreParameterRanges(Gt, rock, ta, 'press_deviation', press_deviation_range);
    [out_resco2]    = exploreParameterRanges(Gt, rock, ta, 'res_co2_sat',     res_co2_sat_range);
    [out_resbri]    = exploreParameterRanges(Gt, rock, ta, 'res_brine_sat',   res_brine_sat_range);
    [out_dismax]    = exploreParameterRanges(Gt, rock, ta, 'dis_max',         diss_max_range);


    % visualization of trends
    % TODO: make following implementation more concise (i.e., put in loop)
    figure; set(gcf,'Position',[100 600 2500 400])
    subplot(1,5,1)
    plot(tgrad_range, out_tgrad.tot_trap_capa_sum, 'bx', ...
        tgrad_range, out_tgrad.strap_mass_co2_sum, 'ro', ...
        tgrad_range, out_tgrad.btrap_mass_co2_res_sum, 'g+', ...
        tgrad_range, out_tgrad.btrap_mass_co2_dis_sum, 'mv');
    legend('total','structural','residual','dissolution')
    xlabel({'Temp Gradient, degrees per kilometer';'(other parameters are default values)'});
    ylabel('Trapping Capacity, Gtons');
    title(name)

    subplot(1,5,2)
    plot(press_deviation_range, out_press.tot_trap_capa_sum, 'bx', ...
        press_deviation_range, out_press.strap_mass_co2_sum, 'ro', ...
        press_deviation_range, out_press.btrap_mass_co2_res_sum, 'g+', ...
        press_deviation_range, out_press.btrap_mass_co2_dis_sum, 'mv');
    legend('total','structural','residual','dissolution')
    xlabel({'Pressure deviation from hydrostatic, percent';'(other parameters are default values)'});
    ylabel('Trapping Capacity, Gtons');
    title(name)

    subplot(1,5,3)
    plot(res_co2_sat_range, out_resco2.tot_trap_capa_sum, 'bx', ...
        res_co2_sat_range, out_resco2.strap_mass_co2_sum, 'ro', ...
        res_co2_sat_range, out_resco2.btrap_mass_co2_res_sum, 'g+', ...
        res_co2_sat_range, out_resco2.btrap_mass_co2_dis_sum, 'mv');
    legend('total','structural','residual','dissolution')
    xlabel({'Residual CO2 saturation';'(other parameters are default values)'});
    ylabel('Trapping Capacity, Gtons');
    title(name)

    subplot(1,5,4)
    plot(res_brine_sat_range, out_resbri.tot_trap_capa_sum, 'bx', ...
        res_brine_sat_range, out_resbri.strap_mass_co2_sum, 'ro', ...
        res_brine_sat_range, out_resbri.btrap_mass_co2_res_sum, 'g+', ...
        res_brine_sat_range, out_resbri.btrap_mass_co2_dis_sum, 'mv');
    legend('total','structural','residual','dissolution')
    xlabel({'Residual Brine saturation';'(other parameters are default values)'});
    ylabel('Trapping Capacity, Gtons');
    title(name)

    subplot(1,5,5)
    plot(diss_max_range, out_dismax.tot_trap_capa_sum, 'bx', ...
        diss_max_range, out_dismax.strap_mass_co2_sum, 'ro', ...
        diss_max_range, out_dismax.btrap_mass_co2_res_sum, 'g+', ...
        diss_max_range, out_dismax.btrap_mass_co2_dis_sum, 'mv');
    legend('total','structural','residual','dissolution')
    xlabel({'Maximum Dissolution';'(other parameters are default values)'});
    ylabel('Trapping Capacity, Gtons');
    title(name)


    % visualization of trends (only structural)
    % TODO: make following implementation more concise (i.e., put in loop)
    figure;  set(gcf,'Position',[100 100 2500 400])
    subplot(1,5,1)
    plot(tgrad_range, out_tgrad.strap_mass_co2_sum, 'ro');
    legend('structural')
    xlabel({'Temp Gradient, degrees per kilometer';'(other parameters are default values)'});
    ylabel('Trapping Capacity, Gtons');
    title(name)

    subplot(1,5,2)
    plot(press_deviation_range, out_press.strap_mass_co2_sum, 'ro');
    legend('structural')
    xlabel({'Pressure deviation from hydrostatic, percent';'(other parameters are default values)'});
    ylabel('Trapping Capacity, Gtons');
    title(name)

    subplot(1,5,3)
    plot(res_co2_sat_range, out_resco2.strap_mass_co2_sum, 'ro');
    legend('structural')
    xlabel({'Residual CO2 saturation';'(other parameters are default values)'});
    ylabel('Trapping Capacity, Gtons');
    title(name)

    subplot(1,5,4)
    plot(res_brine_sat_range, out_resbri.strap_mass_co2_sum, 'ro');
    legend('structural')
    xlabel({'Residual Brine saturation';'(other parameters are default values)'});
    ylabel('Trapping Capacity, Gtons');
    title(name)

    subplot(1,5,5)
    plot(diss_max_range, out_dismax.strap_mass_co2_sum, 'ro');
    legend('structural')
    xlabel({'Maximum Dissolution';'(other parameters are default values)'});
    ylabel('Trapping Capacity, Gtons');
    title(name)


end
                



%% 4. Use method to place wells:
% a) manual -- use interactiveTrapping to study the best region to place
% injection point such that the greastest volume of structurally traps are
% utilized.

if UseInteractiveTrapping

    % color by spill regions. Traps colored grey when colorpath is false,
    % otherwise traps colored red with colorpath is true.
    h = interactiveTrapping(Gt, 'method', 'node', 'light', true, ...
       'spillregions', true, 'colorpath', false); %, 'contour', false) %'injpt', wi_a);

    view(-80,64);
    disp('Showing CO2 Atlas dataset... Close to continue')

    % adjust location of third subplot before saving figure as png
    % note: h.Children(i) where i isn't fixed. Look at h.Children
    disp('Check h.Children')
    pause
    set(h.Children(4), 'position',[0.775 0.09 0.2 0.30])


    % b) use optimization (objective function) to determine well placement of N
    % wells.


    % c) use array of wells

end


%% 5. Given placed wells, compute their injection rates:
% a) compute rate such that an injection period of X years will fill the
% volume of traps along the spill-paths over a period of Y years.
% (user-specified X and Y)




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

% _________________________________________________________________________
% 2) set up schedule (wells, bc, etc.).

% WELLS:
rock2D = rock; % TODO: use por/perm maps
rhoCref = 760 * kilogram / meter ^3; % an (arbitrary) reference density

% TODO: specify coordinate and compute corresponding index of cell (which
% depends on coarsening level)
wellCellIndex = 1554; %[wellCoord_x wellCoord_y wellCoord_z];
[i, j] = ind2sub(Gt.cartDims, wellCellIndex);
wellCoord_x = Gt.cells.centroids(100,1);
wellCoord_y = Gt.cells.centroids(100,2);
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

h = figure;
subplot(1,6,1)
plotCellData(Gt, caprock_temperature)
title('Caprock Temperature'); axis off equal tight
hcb = colorbar;
hcb.Label.String = 'Kelvin'; set(hcb, 'fontSize', 18)

subplot(1,6,2)
plotCellData(Gt, press)
title('Pressure'); axis off equal tight
hcb = colorbar;
hcb.Label.String = 'Pascals'; set(hcb, 'fontSize', 18)

subplot(1,6,3)
plotCellData(Gt, pressDiffFromHydrostatic)
title('Pressure diff. from hydrostatic'); axis off equal tight
hcb = colorbar;
hcb.Label.String = 'Pascals'; set(hcb, 'fontSize', 18)

subplot(1,6,4)
plotCellData(Gt, densityCO2)
title('CO2 density at Caprock'); axis off equal tight
hcb = colorbar;
hcb.Label.String = 'kg/m^3'; set(hcb, 'fontSize', 18)

subplot(1,6,5)
plotGrid(Gt, 'FaceColor', 'none', 'EdgeAlpha', 0.1)
plotCellData(Gt, satCO2, satCO2>(0.1/100)) %satCO2~=0)
title('Saturation of CO2'); axis off equal tight
hcb = colorbar;
hcb.Label.String = 'Pore Volume CO2 Saturation'; set(hcb, 'fontSize', 18)

subplot(1,6,6)
plotGrid(Gt, 'FaceColor', 'none', 'EdgeAlpha', 0.1)
%plotFaces(Gt, bf, 'EdgeColor','k', 'LineWidth',3);
plotCellData(Gt, massCO2/1e9, satCO2>(0.1/100)) % only plot plume that has sat > 0.1 percent
title('Mass of CO2'); axis off equal tight
hcb = colorbar;
hcb.Label.String = 'Mt'; set(hcb, 'fontSize', 18)


%% 7. Use objective function to get optimized injection rates at wells.


