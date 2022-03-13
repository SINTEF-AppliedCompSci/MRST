%% Examples to demonstrate forecast of CO2 leakage and its application
% By assuming a system's flow dynamics are gravity-dominated, we determine
% the amount of CO2 remaining in the domain using spill-point dynamics.

mrstModule add co2lab

%% 1. Synthetic grid
% For simplicity, we first consider a synthetic grid of a sloping and
% perturbed topsurface. Beginning with a domain that is fully saturated
% with CO2, the forecast amount of CO2 remaining is shown to be equivalent
% to the structural trapping capacity (in the event that there is no
% residual trapping). On the other hand, when residual trapping is
% included, the forecast amount remaining is a conservative estimate (since
% we do not forecast the areal sweep of the plume as it migrates), but this
% estimate is shown to be close to the structural and residual trapping
% capacity.

% Getting grid and rock:
Gt = dippedPerturbedGrid();
rock2D.poro = ones(Gt.cells.num,1) * 0.25;
rock2D.perm = ones(Gt.cells.num,1) * 500 * milli*darcy;

%%
% Trap analysis: labeling trap and catchment numbers
figure; hold on
ta = trapAnalysis(Gt, true);
plotGrid(Gt, 'facecolor','none', 'edgealpha',0.1);
colorizeCatchmentRegions(Gt, ta);
num_traps = numel((unique(ta.traps(ta.traps ~= 0))));
for i = 1:num_traps
    trapcells = find(ta.traps==i);
    plotCellData(Gt, ones(Gt.cells.num,1).*i, trapcells, ...
        'facecolor',[0.5 0.5 0.5], 'facealpha',0.5, 'edgecolor','none');
end
top = topCellOfTrap(Gt, ta);
[x, y] = deal(Gt.cells.centroids(top,1), Gt.cells.centroids(top,2));
data = 1:1:numel(top);
data = data';
text(x, y, Gt.cells.z(top)-3, cellstr(num2str(data, '%3.0f')), ...
   'HorizontalAlignment','center', 'FontSize',16, 'FontWeight', 'bold');
axis tight off; view([42,48])

%%
% Fluid model
seainfo = getSeaInfo('NorthSea',760);
%seainfo.res_sat_wat = 0; % @@ can use for testing
%seainfo.res_sat_co2 = 0; % @@ can use for testing
T_ref = computeCaprockTemperature(Gt, seainfo.seafloor_temp, ...
    seainfo.seafloor_depth, seainfo.temp_gradient);
T_ref = T_ref + 273.15; % Kelvin
p     = computeHydrostaticPressure(Gt, seainfo.water_density, 1*atm);
P_ref = mean(p); 
fluid = makeVEFluid(Gt, rock2D, 'sharp interface'                   , ...
       'fixedT'       , T_ref                                       , ...
       'residual'     , [seainfo.res_sat_wat  , seainfo.res_sat_co2], ...
       'wat_rho_ref'  , seainfo.water_density                       , ...
       'co2_rho_ref'  , seainfo.rhoCref                             , ...
       'wat_rho_pvt'  , [4.3e-5/barsa         , P_ref]              , ...
       'co2_rho_pvt'  , [[0.1 400]*mega*Pascal, [4 250]+274]        , ...
       'wat_mu_ref'   , seainfo.water_mu                            , ...
       'co2_mu_pvt'   , [[0.1 400]*mega*Pascal, [4 250]+274]        , ...
       'pvMult_fac'   , 1e-5/barsa                                  , ...
       'pvMult_p_ref' , P_ref                                       , ... 
       'dissolution'  , false                                       , ...
       'surf_topo'    , 'smooth');

%%    
% Creating a state of CO2 saturation:
% We assume the formation is fully saturated with CO2. If there is no
% residual water, then CO2 saturations will be 1, otherwise they are equal
% to 1-sw to account for residual water occupying the pore space.
sG    = ones(Gt.cells.num,1) .* (1 - fluid.res_water);
sGmax = sG;
sF    = ones(Gt.cells.num,1) - sG;

%%
% Using spill-point dynamics to predict future mass remaining/leaked:
mrstVerbose on;
[ will_stay, will_leak ] = massAtInfinity( Gt, rock2D, p, ...
    sG, sGmax, sF, 0, fluid, ta, [], 'plotsOn',true );


%% 2. Injected CO2 scenario into a realistic grid
% Here, we simulate the injection and migration of CO2 into a real
% formation using a vertical-equilibrium model. This example helps to show
% the transition from pressure-driven to gravity-dominated flow, and the
% impact it has on our forecasting algorithm. That is, non-monotonic
% behavior is noticed in the forecast curve immediately following the
% injection period, and several years must be simulated before the forecast
% curve begins to converge. However, once the forecast has converged, there
% is no need to simulate the migration period any further.

% Setting up (cropped) model:
[Gt, rock2D] = getFormationTopGrid('Utsirafm',3);
ind          = (max(Gt.cells.centroids(:,2)) - Gt.cells.centroids(:,2)) < 200 * kilo * meter;
[g, cellmap] = removeCells(Gt.parent, find(~ind)); clear Gt
Gt           = topSurfaceGrid(g);   clear g
rock2D.poro = rock2D.poro(cellmap);
rock2D.perm = rock2D.perm(cellmap); clear cellmap
ta               = trapAnalysis(Gt, true);
seainfo          = getSeaInfo('Utsirafm', 760);
caprock_pressure = (Gt.cells.z * seainfo.water_density * norm(gravity)) ...
                .* (1 + seainfo.press_deviation/100); % @@ should add 1*atm
P_ref = mean(caprock_pressure);
T_ref = computeCaprockTemperature(Gt, seainfo.seafloor_temp, ...
    seainfo.seafloor_depth, seainfo.temp_gradient);
T_ref = T_ref + 273.15; % Kelvin
fluid = makeVEFluid(Gt, rock2D, 'sharp interface'                   , ...
       'fixedT'       , T_ref                                       , ...
       'residual'     , [seainfo.res_sat_wat  , seainfo.res_sat_co2], ...
       'wat_rho_ref'  , seainfo.water_density                       , ...
       'co2_rho_ref'  , seainfo.rhoCref                             , ...
       'wat_rho_pvt'  , [4.3e-5/barsa         , P_ref]              , ...
       'co2_rho_pvt'  , [[0.1 400]*mega*Pascal, [4 250]+274]        , ...
       'wat_mu_ref'   , seainfo.water_mu                            , ...
       'co2_mu_pvt'   , [[0.1 400]*mega*Pascal, [4 250]+274]        , ...
       'pvMult_fac'   , 1e-5/barsa                                  , ...
       'pvMult_p_ref' , P_ref                                       , ... 
       'dissolution'  , false                                       , ...
       'surf_topo'    , 'smooth');               
model = CO2VEBlackOilTypeModel(Gt, rock2D, fluid);

% Defining initial state:
initState.pressure  = computeHydrostaticPressure(Gt, fluid.rhoWS, 1*atm);
initState.s         = repmat([1 0], Gt.cells.num, 1);
initState.sGmax     = initState.s(:,2);

%% Simulate injection and migration
% We simulate 10 years of injection and 300 years of migration using a
% standard VE model.
itime  = 10 * year;
isteps = 10;
mtime  = 300 * year;
msteps = 30;
dTi         = itime / isteps;
dTm         = mtime / msteps;
istepvec    = ones(isteps, 1) * dTi;
mstepvec    = ones(msteps, 1) * dTm;
schedule.step.val       = [istepvec; mstepvec];
schedule.step.control   = [ones(isteps, 1); ones(msteps, 1) * 2];

% Setting up injection well:
wcinx = findEnclosingCell(Gt, [4.7777e5, 6.8154e6]);
W = addWell([], Gt.parent, rock2D, wcinx, ...
                'name',     ['Winj' num2str(wcinx)],  ...
                'Type',     'rate', ...
                'Val',      0.3536, ...
                'comp_i',   [0,1], ...
                'Radius',   0.3 );
W = convertwellsVE(W, Gt.parent, Gt, rock2D, 'ip_tpf');
W_shut = W;
W_shut.val = sqrt(eps); 
schedule.control(1).W = W;
schedule.control(2).W = W_shut;

% Setting up boundary conditions:
bfaces  = find(any(Gt.faces.neighbors==0,2));
bdryVal = Gt.faces.z(bfaces) * fluid.rhoWS * norm(gravity) + 1*atm;
schedule.control(1).bc = addBC( [], bfaces, 'pressure', bdryVal, 'sat', [1 0] );
schedule.control(2).bc = addBC( [], bfaces, 'pressure', bdryVal, 'sat', [1 0] );

% Simulating injection and migration:
mrstVerbose off;
[wellSols, states] = simulateScheduleAD(initState, model, schedule);

%% Forecast
% For each time step in the simulation, we extract the aquifer state and
% use it to make a spill-point forecast of the mass that will be trapped at
% time infinity. We then plot the resulting curve as a blue line on top of
% the trapping inventory from the vertical-equilibrium simulation.
Ma = zeros(numel(states),1);
for i = 1:numel(states)
    Ma(i) = massAtInfinity( model.G, model.rock, states{i}.pressure, ...
        states{i}.s(:,2), states{i}.sGmax, states{i}.s(:,1), 0, ...
        model.fluid, ta, [] ) / 1e12; % Gt
end

% Plotting Ma on top of trapping inventory
reports = makeReports(Gt, [{initState}; states], rock2D, fluid, ...
                    schedule, [fluid.res_water, fluid.res_gas], ta, []);
h = figure;
plot(1);
ax = get(h,'currentaxes');
plotTrappingDistribution(ax, reports, 'legend_location', 'southeast');
fsize = 10;
set(get(gca, 'xlabel'), 'fontsize', fsize)
set(get(gca, 'ylabel'), 'fontsize', fsize)
set(gca,'fontsize', fsize);
hold on;
plot([0; cumsum(convertTo(schedule.step.val,year))], [0; Ma*1e3], 'b','LineWidth',3)

%% Show non-monotonicity
% Adding stars to forecast curve corresponding to the CO2 sat. snapshots
[x,y] = deal([0; cumsum(convertTo(schedule.step.val,year))], [0; Ma*1e3]);
% x is in years, y is in Mt
for i=1:5
    plot(x(10+i), y(10+i), '*', 'MarkerSize',10, 'LineWidth',2, 'Color','k')
end

% Zooming into curve to see non-monotonicity
xlim([0 100])
ylim([83 84.738])

%% Plot snapshots of CO2 saturation
sts = [{initState}; states];
tol = 0.005;
time_yr = [0; convertTo(cumsum(schedule.step.val),year)];
for np = 1:5
    
    tinx = 10 + np; % plot "sts" corresponding to years 10, 20, 30, 40, 50

    figure; hold on
    plotGrid(Gt, 'facecolor','none', 'edgealpha',0.1)
    colorizeCatchmentRegions(Gt, ta);
    plotCellData(Gt, sts{tinx}.s(:,2), sts{tinx}.s(:,2) > tol, 'edgecolor','none')

    cmax = max(sts{tinx}.s(:,2));
    set(findobj(gcf,'type','axes'),'clim',[0, cmax])
    colormap(gca,parula)
    daspect([1 1 0.02])
    axis tight off

    % add rivers between traps
    rivers = ta.cell_lines;
    for tr = rivers
        for r = tr{:}
            drawCellConnections3D(Gt, r{:}, 'color', 'k', 'lineWidth', 2);
        end
    end

    % add traps (or just trap outlines)
    plotCellData(Gt, ones(Gt.cells.num,1), ta.traps ~= 0, ...
        'facecolor','b','facealpha',0.2, 'edgecolor','none')
    for i = 1:numel(ta.trap_z)
        plotFaces(Gt, boundaryFaces(Gt, ta.trap_regions == i), 'EdgeColor','k','LineWidth',1)
        plotFaces(Gt, boundaryFaces(Gt, ta.traps == i), 'EdgeColor','b','LineWidth',2)
    end

    % add well location
    plotWell(Gt.parent, W, 'color','k')
    
    
    title(['Year ',num2str(time_yr(tinx))])
    
    % getting reports for catchment regions surrounding the well
    fprintf('\n-----------\n')
    fprintf(  'By Year %d:', time_yr(tinx))
    fprintf('\n-----------\n')
    massAtInfinity( Gt, rock2D, sts{tinx}.pressure, sts{tinx}.s(:,2), ...
        sts{tinx}.sGmax, sts{tinx}.s(:,1), 0, model.fluid, ta, [], ...
        'plotsOn',false , ...
        'report_trap_regions',[37, 45, 48, 53, 52], ...
        'tot_inj',sum(reports(tinx).masses));

end

%%
% <html>
% <p><font size="-1">
% Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.
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
