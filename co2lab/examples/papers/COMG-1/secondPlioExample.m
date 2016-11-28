%% Inspect the Pliocenesand Formation
% This formation has a very low percentage (0.02%) of trapping compared to
% the overall volume of the whole model

mrstModule add coarsegrid deckformat mex ad-core ad-props
grdecl = getAtlasGrid('Pliocenesand');
G      = processGRDECL(grdecl{1});
G      = computeGeometry(G(1));
Gt     = topSurfaceGrid(G);

%%
% Next, we will perform a VE simulation. This could, of course, have been
% launched from inside the interactive viewer, but to make the example as
% reproducible as possible, we launch it manually from the outside.

%%
% Rerun for a longer time
petrodata.avgperm = 1.2*darcy;
petrodata.avgporo = 0.25;

%%
% Adding depth to the plio example
% It is to shallow for a real storage site
depth=1200;

% cut grid to avoid calculation on not relevant domain
wpos = Gt.parent.cells.centroids(5280, 1:2);
wpos(:, 1) = 4.85e5;
G = Gt.parent;
rm_cells = abs(Gt.cells.centroids(:, 2) - wpos(:, 2))>2.5e4;
G = removeCells(G, rm_cells);
G.nodes.coords(:, 3) = G.nodes.coords(:, 3) + depth;
G = computeGeometry(G);
Gt = topSurfaceGrid(G);

% set the permeability and porosity
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


Tm  = 3000*year;
dTm = 25*year/2; 

mstep = linspace(0.5*year, dTm, 5)';
mstep = [mstep; ones(floor((Tm-sum(mstep))/dTm),1)*dTm];
mstep = [mstep; Tm-sum(mstep)];

%%
amount = 5; %Mt/year

% convert mass flux at surface to surface reference volume/"surface volume
rhoc = 760; rhow = 1100;
rate = amount * 1e9 * kilogram ./ (year * rhoc * kilogram *meter^3);

% Specify residual saturations
res_water = 0.11;
res_gas = 0.21;

% find well position
dist = sqrt(sum(bsxfun(@minus, Gt.cells.centroids(:, 1:2), wpos).^2, 2));
[dd, cellnum] = min(dist);
[ix, iy] = ind2sub(Gt.cartDims, Gt.cells.indexMap(cellnum));
wellIx = double([ix iy]);

% make well
assert(G.cartDims(3)==1)
W = createSampleWell([],Gt.parent, rock, cellnum ,     ...
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
k = 1;
for dis_model = {'none', 'instant', 'rate'};
   % Set up fluid model
   p_range = [0.1, 400] * mega * Pascal; % CO2 default pressure range
   t_range = [  4, 250] + 274;           % CO2 default temperature range
   cw      = 4.3e-5 / barsa; % water linear compressibility
   temp_grad = 30 / (kilo*meter); % temperature gradient
   surf_temp = 12; % surface temperature in Celsius
   temperature = Gt.cells.z .* temp_grad + (274 + surf_temp); % fixed temperature

   fluid = makeVEFluid(Gt, rock2D, 'sharp interface', ...
                       'residual', [res_water, res_gas], ...
                       'co2_rho_ref', rhoc , ...
                       'wat_rho_ref', rhow , ...
                       'co2_mu_ref', 6e-5 , ...
                       'wat_mu_ref', 8e-4 , ...
                       'co2_rho_pvt', [p_range, t_range], ...
                       'wat_rho_pvt', [cw, 100 * barsa], ...
                       'fixedT', temperature, ...
                       'dissolution', ~strcmpi(dis_model{:}, 'none'), ...
                       'dis_rate', strcmpi(dis_model{:}, 'rate') * 5e-11, ...
                       'dis_max' , 0.03);

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
   ensure_path_exists('data/');
   save(['data/secondPlioExample_',num2str(depth),'_',num2str(k),'.mat'], ...
        't2','states','wellSols','schedule', 'Gt', 'fluid', 'rock2D', 'report');
   k = k+1;
end
