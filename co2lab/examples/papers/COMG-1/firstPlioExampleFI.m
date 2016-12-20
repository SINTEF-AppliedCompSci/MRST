%% Pliocenesand Formation: fully-implicit simulation
% We run the first Pliocenesand example using a fully implicit simulator to
% show the improved computational efficiency compared with the sequentially
% implicit formulation

mrstModule add coarsegrid deckformat mex ad-core ad-props
grdecl = getAtlasGrid('Pliocenesand');
G      = processGRDECL(grdecl{1});
depth=1200;
G.nodes.coords(:, 3) = G.nodes.coords(:, 3) + depth;
G = computeGeometry(G);
Gt = topSurfaceGrid(G);

% set the permeability and porosity
petrodata.avgperm = 1.2*darcy;
petrodata.avgporo = 0.25;
rock.poro = repmat(petrodata.avgporo, Gt.cells.num, 1);
rock.perm = repmat(petrodata.avgperm, Gt.cells.num, 1);
rock2D    = averageRock(rock, Gt);

%% Find pressure boundary
% Setting boundary conditions is unfortunately a manual process and may
% require some fiddling with indices, as shown in the code below. Here, we
% need to find all outer vertical faces
i = any(Gt.faces.neighbors==0, 2);  % find all outer faces
I = i(Gt.cells.faces(:,1));         % vector of all faces of all cells, true if outer
j = false(6,1);                     % mask, cells can at most have 6 faces,

j(1)=true;
%   extract east, west, north, south
J = j(Gt.cells.faces(:,2));         % vector of faces per cell, true if E,W,N,S
bcIxVE = Gt.cells.faces(I & J, 1);

%% Set time and fluid parameters

% Fluid data are taken from paper SPE 134891
gravity on

Ti  =   50*year;

dTi =  2*year;
istep = linspace(0.1*year, dTi, 10)';
istep = [istep; ones(floor((Ti-sum(istep))/dTi), 1)*dTi];
istep = [istep; Ti-sum(istep)];


Tm  = 1450*year;
dTm = 10*year; 

mstep = linspace(0.5*year, dTm, 5)';
mstep = [mstep; ones(floor((Tm-sum(mstep))/dTm),1)*dTm];
mstep = [mstep; Tm-sum(mstep)];

%%
amount = 10; %Mt/year

% convert mass flux at surface to surface reference volume/"surface volume
rhoc = 686.54; rhow = 975.86;
rate = amount * 1e9 * kilogram ./ (year * rhoc * kilogram *meter^3);

% Specify residual saturations
res_water = 0.1;
res_gas = 0.2;

% find well position
coord = [464328, 6646937]; 
dist2 = sum(bsxfun(@minus, Gt.cells.centroids, coord).^2, 2);
[~, cellnum] = min(dist2);
[ix, iy] = ind2sub(Gt.cartDims, Gt.cells.indexMap(cellnum));
wellIx = double([ix iy]);

% make well
assert(G.cartDims(3)==1)
W = createSampleWell([],Gt.parent, rock, cellnum,     ...
                        'Type', 'rate', 'Val', rate, ...
                        'Radius', 0.125, 'Name', 'I','Comp_i',[0 1]);
W = convertwellsVE(W, Gt.parent, Gt, rock2D, 'ip_tpf');
W_shut = W;
W_shut.val = 0;

% specify boundary conditions
bc = addBC([], bcIxVE, 'pressure', Gt.faces.z(bcIxVE) * rhow * norm(gravity));
bc.sat = [ones(numel(bc.face), 1), zeros(numel(bc.face), 1)];

% make schedule
schedule.control = [struct('W', W, 'bc', bc), struct('W', W_shut, 'bc', bc)];
schedule.step = struct('control', [ones(size(istep));ones(size(mstep))*2], ...
                       'val', [istep; mstep]);

%% start loop over cases
% Set up fluid model
p_range = [0.1, 400] * mega * Pascal; % CO2 default pressure range
t_range = [  4, 250] + 274;           % CO2 default temperature range
cw      = 1e-8/barsa; %4.3e-5 / barsa; % water linear compressibility
temp_grad = 30 / (kilo*meter); % temperature gradient
surf_temp = 12; % surface temperature in Celsius
temperature = Gt.cells.z .* temp_grad + (274 + surf_temp); % fixed temperature

fluid = makeVEFluid(Gt, rock2D, 'sharp interface', ...
                    'residual', [res_water, res_gas], ...
                    'co2_rho_ref', rhoc , ...
                    'wat_rho_ref', rhow , ...
                    'co2_mu_ref', 0.056641*centi*poise, ...
                    'wat_mu_ref', 0.30860*centi*poise , ...
                    'co2_rho_pvt', [cw, 100*barsa], ...
                    'wat_rho_pvt', [cw, 100 * barsa], ...
                    'fixedT', temperature, ...
                    'dissolution', false);
                
fluid = addVERelperm(fluid       , Gt          , ...
                     'res_water' , res_water , ...
                     'res_gas'   , res_gas , ...
                     'krg',       0.2142 , ...
                     'krw',       0.85 );
                 
% Set up initial conditions
nc = Gt.cells.num;
initState = struct('pressure', Gt.cells.z(:)*norm(gravity)*fluid.rhoWS, ...
                   's'       , [ones(nc, 1), zeros(nc, 1)], ...
                   'sGmax'   , zeros(nc, 1), ...
                   'rs'      , zeros(G.cells.num, 1));

% Setup model and run simulation
model = CO2VEBlackOilTypeModel(Gt, rock2D, fluid);
t2 = tic;
[wellSols, states, report] = ...
    simulateScheduleAD(initState, model, schedule, ...
                       'NonLinearSolver', NonLinearSolver('useRelaxation', true));
t2 = toc(t2);

states = {initState, states{:}}';%#ok  % Include initial state
k=1;
ensure_path_exists('data/');
save(['data/secondPlioExample_',num2str(depth),'_',num2str(k),'.mat'], ...
    't2','states','wellSols','schedule', 'Gt', 'fluid', 'rock2D', 'report');

