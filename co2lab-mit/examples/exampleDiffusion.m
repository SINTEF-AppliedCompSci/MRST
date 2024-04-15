%% Perform CO2 injection (blackoil module)
% MRST Reference: Lie, K.-A., 2019. An Introduction to Reservoir Simulation 
% using MATLAB: <https://doi.org/10.1017/9781108591416 https://doi.org/10.1017/9781108591416> 
% (free PDF)
% 
% See also:  <https://www.sintef.no/projectweb/mrst/ https://www.sintef.no/projectweb/mrst/>
% 
% Download MRST:  <https://www.sintef.no/projectweb/mrst/download/ https://www.sintef.no/projectweb/mrst/download/>
% 
% Example prepared by Lluis Salo-Salgado (<mailto:lsalo@mit.edu lsalo@mit.edu>)
% 
% How to run this script: 
% 
% 1. Download MRST (see link above) 
% 
% 2. Open MATLAB 
% 
% 3. Run the startup.m file within the MRST folder 
% 
% 4. You can now run the script. Note that, the first time you run the model, 
% you'll have to accept the prompts to download and setup AMGCL and boost. MRST 
% takes care of this automatically so you just need to click "OK". This is required 
% to run the fine mesh and large models in general, otherwise Matlab is too slow.
% 
% What does this code do? Run a CO2 injection simulation with a quasi-2D mesh. 
% We initialize the model with 3 different rock/saturation regions (reservoir, 
% caprock/seal and fault) and full water saturation. We then inject CO2 for a 
% while and let the model run. We consider capillary pressure and dissolution 
% of CO2 in the water. This is a blackoil-type model, so the "gas" phase is the 
% CO2-rich phase, and the "oil" phase is given the properties of water/brine.

% Cleanup and load modules. You will not need ls-proj.
clear, close all
mrstModule add upr ad-props ad-blackoil deckformat ad-core mrst-gui linearsolvers  
mrstVerbose on

%% Output directory and options

% Options
plotFigs = true;
mesh    = 'coarse';                           % 'coarse', 'medium', or 'fine'
wellno  = 1;                                  % 1 or 2 injectors (TBD for 2)
D       = 2e-11;                              % pseudo-molecular diffusivity (m^2/s)
rate    = 8;                                  % max inj rate (mL/min, surface conditions) 
folderName = ['example_mesh_', mesh, '_rate_', num2str(rate), '_depth_1km_D_' num2str(D)];

% Plots
fig3D = @() figure('Position', [0, 0, 1300, 650]);
alpha = 0.6;
cmap  = jet*alpha + (1-alpha);
setAxProps = @(ax) set(ax, 'View'              , [65, 20]         , ...
                           'PlotBoxAspectRatio', [4.40, 1.86, 1.00], ...
                           'Projection'        , 'Perspective'     , ...
                           'Box'               , 'on'              , ...
                           'ColorMap'          , cmap              );


%% 1. Mesh
% Load simpleExtrudedMesh
G_dat = simpleExtrudedFluidFlowerMesh(mesh);
sealID = [1 3 5 8];
faultID = 6;
unitGroups = {5, [4 7], [1 8], [2 9], 3, 6};    % bot to top, fault at the end
G = G_dat.G;

% Plot full grid
if plotFigs
    fig3D();
    %plotCellData(G, G_dat.p)
    Ncompart = max(G_dat.p);
    cmap = copper(Ncompart);
    colr = [1, 7, 4, 8, 2, 5];
    for n=1:numel(unitGroups)
        plotCellData(G, colr(n)*ones(sum(ismember(G_dat.p, unitGroups{n})), 1), ...
                     ismember(G_dat.p, unitGroups{n}), 'edgealpha', 0.2)
    end
    plotGrid(G, find(G_dat.wellNo==1), 'facecolor', 'r')
    outlineCoarseGrid(G, G_dat.compartID,'EdgeColor','w','LineWidth',2);
    setAxProps(gca), %camlight();
    colormap(copper); c = colorbar; set(c, 'YTick', sort(colr));
    axis equal on
    %ylim([0 1]), zlim([0 1])
end

% Remove cells of seal, since we are interested in convective mixing in reservoir
%  and migration along potential fault conduit
idc = 1:G.cells.num;
G.faces = rmfield(G.faces, 'tag');
[G, cellmap] = removeCells(G, ismember(G_dat.compartID, sealID));


%% 2. Rock
% Here we define the porosity and permeability [L^2] of each grid cell. Layers 
% follow unitGroups (bot to top, left to right for faulted, fault at the end)
% find fault cell ids
cid = 1:G_dat.G.cells.num;
fid = cid(G_dat.compartID == faultID);
fid = ismember(cellmap, fid);
%plotGrid(G,'facecolor','none'); plotGrid(G, fid); view([90 0]);

% populate poro perm
rock.poro = 0.3*ones(G.cells.num, 1);
rock.poro(fid) = 0.1;
rock.perm = 1000*ones(G.cells.num, 1);    % isotropic
rock.perm(fid) = 10;
rock.perm = rock.perm*(milli*darcy);

% Plot permeability
if plotFigs
    fig3D(); plotGrid(G_dat.G, ismember(G_dat.p, sealID), 'facecolor', 'none', 'edgealpha', 0.1);
    plotCellData(G, log10(rock.perm/(milli*darcy)), 'edgealpha', 0.2)
    setAxProps(gca), colormap(copper), 
    if strcmp(mesh, 'fine')
        plotGrid(G, 27110, 'facecolor', 'r', 'edgecolor', 'none')
    end
    c = colorbar; c.Label.Interpreter = 'latex'; c.Label.FontSize = 11;
    c.Label.String = '$\log_{10}(k $ [mD])';
    axis equal off; view ([90 0])
end


%% 3. Fluids
% Here we define the fluid properties. This is a blackoil model, so we need 
% a fluid object with the corresponding formation volume factors (FVF) and solution 
% gas oil ratio (Rs). We obtained the PVT properties using the pvtBrineWithCO2BlackOil 
% function (co2-brine_props) for a saline aquifer at ~1km depth.
% 
% In addition, we need the relative permeability, capillary pressure, surface 
% density of the fluids and rock compressibility. The fluid object containing 
% these properties (variable fluid) will have 3 cells for those properties changing 
% from one model region to another (reservoir, caprock and falut).

% Define saturation regions
rock.regions.saturation = ones(G.cells.num, 1);        % reservoir units
%rock.regions.saturation(ismember(G_dat.compartID, sealID)) = 2; % seal no longer present
rock.regions.saturation(fid) = 3;

% Define 1 region for compressibility
rock.regions.rocknum = ones(G.cells.num,1); 

% Load fluid deck from .DATA (ECLIPSE-type) input file and initialize fluid
fluid_path = fullfile(getDatasetPath('co2labmit'), 'fluid_props');
fn  = fullfile(fluid_path, 'example_co2brine_1kmDepth_3regions.DATA');
deck = convertDeckUnits(readEclipseDeck(fn));
deck.REGIONS.ROCKNUM = rock.regions.rocknum;
fluid = initDeckADIFluid(deck);                 % this is the fluid object

% Plot fluid density and viscosity
if plotFigs
    np          = 50;
    p_val       = linspace(80,160,np)'*barsa;    % MRST always SI units (Pa)
    rho_co2     = fluid.rhoGS*fluid.bG(p_val);   % Dry gas (no water in gas phase)
    mu_co2      = fluid.muG(p_val);
    rss_val     = fluid.rsSat(p_val);            % Live oil (aqueous phase with dissolved CO2)
    rho_b_sat   = fluid.bO(p_val,rss_val,true(np,1)) .* ...
                          (rss_val.*fluid.rhoGS + fluid.rhoOS);
    rho_b       = fluid.rhoOS*fluid.bO(p_val,zeros(np,1),false(np,1));
    mu_b_sat    = fluid.muO(p_val,rss_val,true(np,1));
    mu_b        = fluid.muO(p_val,zeros(np,1),false(np,1));

    % CO2
    figure(11)
    tiledlayout(1,2,"TileSpacing","compact","Padding","compact")
    nexttile(1)
    hold on
    plot(p_val/barsa, rho_co2, '-r', 'linewidth', 1.5); 
    plot(p_val/barsa, rho_b, '-c', 'linewidth', 1.5); 
    plot(p_val/barsa, rho_b_sat, '-b', 'linewidth', 1.5); 
    hold off
    grid on
    xlabel('p [bar]', 'fontsize', 12)
    ylabel('\rho [kg/m^3]', 'fontsize', 12)
    xlim([80 160])
    ylim([0 1100]),
    legend('Gas', 'water', 'water+CO2 (sat.)', 'location', 'best')
    title('Density', 'fontsize', 14)
    nexttile(2)
    hold on
    plot(p_val/barsa, mu_co2*1e3, '-r', 'linewidth', 1.5); 
    plot(p_val/barsa, mu_b*1e3, '-c', 'linewidth', 1.5); 
    plot(p_val/barsa, mu_b_sat*1e3, '-b', 'linewidth', 1.5); 
    hold off
    grid on
    xlabel('p [bar]', 'fontsize', 12)
    ylabel('\mu [cP]', 'fontsize', 12)
    xlim([80 160])
    ylim([0.01 1])
    set(gca,'Yscale','log')
    title('Viscosity', 'fontsize', 14)
end


%% 4. Initialize (pressure, saturation, rs, rv)
% We compute the hydrostatic distribution numerically, relative to fixed datum 
% point p(z0) = p_r. See MRST book, page 206. Saturation is initialized assuming 
% the full model is water-saturated rs
gravity reset on
g = norm(gravity);
rho_wr = fluid.rhoOS*kilogram/meter^3;
water_column = 1000;                   % 1km depth
p_r = 1*barsa + g*rho_wr*water_column; % p at shallowest z
[z_0, z_max] = deal(min(G.cells.centroids(:,3)), max(G.cells.centroids(:,3)));
equil  = ode23(@(z,p) g .* fluid.bO(p,0,false)*fluid.rhoOS, [z_0, z_max], p_r);
p0 = reshape(deval(equil, G.cells.centroids(:,3)), [], 1);  clear equil
s0  = repmat([1, 0], [G.cells.num, 1]);  % s: fully saturated in oil --> 
                                         %    [0 1 0] if 'WOG'; [1 0] if 'OG'
rs0 = zeros(G.cells.num, 1);             % no dissolved gas at the beginning
rv0 = 0;                                 % dry gas
state0 = struct('s', s0, 'rs', rs0, 'rv', rv0, 'pressure', p0);

% plot
if plotFigs
    fig3D(); plotCellData(G, p0/barsa, 'edgealpha', 0.2)
    setAxProps(gca), colormap(jet), c = colorbar;
    c.Label.Interpreter = 'latex'; c.Label.FontSize = 11;
    c.Label.String = '$p_0 $ [bar]';
    axis equal %off
end


%% 5. Wells
% All times are defined in sec in MRST, so we convert them to sec.
t = [60*minute 1*day 30*day];
reportTimes = [(12:12:t(1)/minute)*minute, ... % ramp up
               (2:1:24)*hour, ...              % injection
               (1440+5:5:1465)*minute, ...     % ramp down
               ([25 26 28 32 36 40 48 60 72 96 120])*hour, ...
                (6:30)*day];
injrate = rate*(milli*litre)/(minute*wellno);    % volumetric rate in [Sm^3/s]
%injvol = injrate*t(2);
%rhoInj = fluid.rhoGS;
%mrate = injvol*rhoInj/(t(2)*wellno);            

%wellInx = find(G_dat.wellNo==1);                 % Well cell id
%wellInx = find(cellmap == wellInx);
dist = pdist2(G.cells.centroids, [0.005, 0.903, 1000.58]);
[~, wellInx] = min(dist,[],1);
W = addWell([ ], G, rock, wellInx, 'Name', 'I1', 'Dir', 'z', ...
            'Type', 'rate', 'Val', injrate, 'compi', [0, 1], ...    % order 'OG'
            'refDepth', G.cells.centroids(wellInx, G.griddim), ...
            'Radius', 1e-3);                                
timesteps = [reportTimes(1) diff(reportTimes)];
assert(sum(timesteps)==t(end), 'sum of timesteps must equal simTime')


%% 6. Model
% Here we put together the blackoil model, and define the acceleration and solver 
% parameters.
model = GenericBlackOilModel(G, rock, fluid, 'disgas', true, 'water', false);
                         
% Acceleration and solver parameters
model.AutoDiffBackend = DiagonalAutoDiffBackend('useMex', true, ...
                                                'deferredAssembly', true);
%model.toleranceCNV    = 2e-3;
%model.dsMaxAbs        = 0.05;
model.minimumPressure = min(state0.pressure);
model = model.validateModel();    

% Diffusion
if D > 0
    diffFlux = CO2TotalFluxWithDiffusion(model);
    diffFlux.componentDiffusion = [0 D];  % opt.D = diffusion coeff.
    diffFlux.faceAverage = true;
    model.FlowDiscretization.ComponentTotalFlux = diffFlux;
end

% nonlinear solver
nls = getNonLinearSolver(model, 'TimestepStrategy', 'iteration', ...
                         'useCPR', true);
nls.useLinesearch = true;
nls.maxIterations = 10;      
nls.maxTimestepCuts = 12; 
nls.acceptanceFactor = 2;


%% BCs
% We impose no-flow everywhere, but we recreate an open aquifer system by multiplying 
% the pore volumes of external cells by a very large number.

% Find external faces
L = max(G.faces.centroids(:,2));
B = max(G.faces.centroids(:,3));
f = any([G.faces.centroids(:,2) == 0, ...
         G.faces.centroids(:,2) > L-1e-3], 2);
%fig3D(); plotFaces(G); plotFaces(G, f, 'edgecolor', 'r', 'facecolor', 'none')
%setAxProps(gca), colormap(jet), axis equal off

% Find external cells
cellsext = unique(reshape(G.faces.neighbors(f, :), [], 1));
cellsext(cellsext==0) = [];
model.operators.pv(cellsext) = model.operators.pv(cellsext)*10^5;
bc = [];

% plot
if plotFigs
    fig3D(); plotGrid(G); hold on, plotGrid(G, cellsext, 'facecolor', 'r')
    setAxProps(gca), axis equal off, view([90 0])
end


%% Schedule
% Here we set up the injection schedule, with proper controls so that we can 
% have changing injection rates.

schedule_inj = simpleSchedule(timesteps, 'W', W, 'bc', bc);      
n_ramp = 5;
v = injrate;
injrates = [0.01*v 0.1*v 0.2*v 0.5*v v ...      % we ramp up in n_ramp steps
            0.8*v 0.5*v 0.1*v 0.01*v 0];
tmp = cell(numel(injrates), 1);                 % create 2 schedules
schedule = struct('step', schedule_inj.step);   % timesteps and wells for each timestep
schedule.control = struct('W', tmp, 'bc', tmp, 'src', tmp); 

% Update injection rates
for n=1:numel(injrates)
    schedule.control(n).W = W;                   % make a copy of well params
    schedule.control(n).W.val = injrates(n);     % update injrate
    schedule.control(n).bc = bc;                 % same BCs
end

% Finally, indicate the index of the well to be used in each timestep
idStep = find(cumsum(schedule.step.val) < t(1));
schedule.step.control(idStep) = 1:max(idStep);
schedule.step.control(idStep(end)+1:end) = max(idStep)+1;
idStep2 = find(cumsum(schedule.step.val) > t(2), 1);
schedule.step.control(idStep2:idStep2+(n_ramp -1)) = (1:n_ramp) + n_ramp;
schedule.step.control(idStep2+n_ramp:end) = 2*n_ramp;


%% Simulation
% For each timestep, MRST will save three files in the folder indicated in outputDir: 
% 
% - report: contains info on the number of iterations and convergence, etc 
% 
% - states: solution structure with p, saturation, etc 
% 
% - wellSols: well parameters such as bhp and such 
% 
if strcmp(mesh, 'fine') || strcmp(mesh, 'medium')
    N = 1;  
    maxNumCompThreads(N);   
end
if isfield(nls.LinearSolver, 'amgcl_setup.nthreads')
    nls.LinearSolver.amgcl_setup.nthreads = N;   % Specify threads manually
end
problem = packSimulationProblem(state0, model, schedule, folderName, ...
                                'Name', folderName, ...
                                'NonLinearSolver', nls);
[ok, status] = simulatePackedProblem(problem);
[wellSols, states, report] = getPackedSimulatorOutput(problem);

%% Results
% Compute quantitites
model = model.validateModel();
for n=1:numel(states)
    states{n}.dp = states{n}.pressure - state0.pressure;  
    pc = model.getProp(states{n}, 'CapillaryPressure');
    states{n}.FlowProps.CapillaryPressure = pc{1,2};
    %states{n}.FlowProps.RelativePermeability = model.getProp(states{n}, ...
    %                                                 'RelativePermeability');
end
states{1}.reg = model.rock.regions.saturation;

% Basic overview
fig3D(); plotToolbar(G, states, 'edgealpha', 0.2); 
setAxProps(gca), colormap(turbo), c = colorbar; %clim([0 40000])
axis equal off
view([90 0])
%set(gca, 'ColorScale', 'log')
%caxis([1e-8 1])

% Flux velocities
f_wat_vol = cell2mat(cellfun(@(x) x.flux(:,1), states, 'uniformoutput', false)'); % volumetric (m3/s)
f_wat = f_wat_vol ./ G.faces.areas;
f_wat_max = max(f_wat);
idv     = G.faces.normals(:,3) == 0;
idh     = all([G.faces.normals(:,1)==0, G.faces.normals(:,2)==0], 2);
f_wat_h = f_wat_vol(idv,:) ./ G.faces.areas(idv);
f_wat_v = f_wat_vol(idh,:) ./ G.faces.areas(idh);
f_wat_h_max = max(f_wat_h);
f_wat_v_max = max(f_wat_v);