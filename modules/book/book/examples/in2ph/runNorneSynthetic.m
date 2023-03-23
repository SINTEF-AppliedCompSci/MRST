%% Two-Phase Incompressible Simulation on the Norne Model
% In this example, we will consider the grid and the petrophysical
% properties from the Norne simulation model to demonstrate some of the
% complexity one will see in a real simulation model. We set up a two-phase
% simulation that uses synthetic fluid model, well placement, and
% simulation schedule.
%
% Norne is an oil and gas field lies located in the Norwegian Sea. The
% reservoir is found in Jurrasic sandstone at a depth of 2500 meter below
% sea level. Operator Statoil and partners (ENI and Petoro) have agreed
% with NTNU to release large amounts of subsurface data from the Norne
% field for research and education purposes.  The
% <http://www.ipt.ntnu.no/~norne/wiki/doku.php Norne Benchmark> datasets
% are hosted and supported by the Center for Integrated Operations in the
% Petroleum Industry (IO Center) at NTNU. Recently, the
% <http://www.opm-project.org OPM Initiative> released the full simulation
% model as an open data set on <https://github.com/OPM/opm-data GitHub>.

mrstModule add deckformat incomp
gravity reset on

%% Read and process the model
% As the Norne dataset is available from the OPM project's public GitHub
% repository, we can download a suitable subset of the simulation model and
% process that subset. See the <matlab:edit('showNorne.m')
% showNorne> script for a more detailed walkthrough of the model.

if ~ (makeNorneSubsetAvailable() && makeNorneGRDECL())
   error('Unable to obtain simulation model subset');
end
grdecl = fullfile(getDatasetPath('norne'), 'NORNE.GRDECL');
grdecl = readGRDECL(grdecl);
grdecl = convertInputUnits(grdecl, getUnitSystem('METRIC'));

%% Extract the active part of the model
% To get the whole grid, we need to reprocess the input data. This gives a
% grid with two components: the main reservoir and a detached stack of
% twelve cells, which we ignore henceforth.
G = processGRDECL(grdecl);
G = computeGeometry(G(1));

%% Set view parameters and plot grid
myview = struct(...
    'view',     [-110,50],   ...  % view angle
    'zoom',     1.8,         ...  % zoom
    'dataasp',  [1 1 .1],    ...  % data aspect ratio
    'dolly',    [0 .1 0],    ...  % camdolly
    'pargs',    {{'EdgeAlpha'; 0.1; 'EdgeColor'; 'k'}} ...
    );
clf
plotReservoirModel(G, [], [], myview);

%% Outline petrophysical properties
% The petrophysical properties are included in the simulation model subset
% represented by function |makeNorneSubsetAvailable|.  Consequently, the
% necessary data has been read into the |grdecl| structure.  We can then
% extract the petrophysical properties. Notice that here, the rock
rock = grdecl2Rock(grdecl, G.cells.indexMap);

%% Porosity and net-to gross
% Show the porosities mapped onto the structural grid
figure,
plotReservoirModel(G, rock.poro, [], myview);
colorbarHist(rock.poro, [.05 .35],'South',100);

%%
% show also the net-to-gross
figure
plotReservoirModel(G, rock.ntg, [], myview);
colorbarHist(rock.ntg,[.07 1],'South',100);


%% Permeability
% The permeability is generally a tridiagonal tensor K = diag(Kx, Ky, Kz).
% For the Norne model, the data file only specifies Kx, and then various
% manipulations are done to get the correct vertical permeability.
figure, 
prms      = myview; 
prms.cs   = [1 10 100 1000 10000];
prms.unit = milli*darcy;
prms.cb   = 'South';
plotReservoirModel(G, log10(rock.perm(:,1)), [], prms);

%%
figure, 
prms.cs   = [0.1 1 10 100 1000 3000];
plotReservoirModel(G, log10(rock.perm(:,3)), [], prms);
%caxis(log10([.1 2500]*milli*darcy));

%% Multipliers
% The model has a number of multipliers that reduce the transmissibility in
% the z-direction
figure
prms = myview;
prms.pargs = {'FaceColor','none','EdgeAlpha',.1,'EdgeColor','k'};
plotReservoirModel(G, [], [], prms);
mz = grdecl.MULTZ(G.cells.indexMap);
plotCellData(G,mz,mz<1,myview.pargs{:});
mz(mz==1) = nan;
colorbarHist(mz, [-.05 .6], 'South', 100);

%% Introduce wells
% The reservoir is produced using a set production wells controlled by
% bottom-hole pressure and rate-controlled injectors. Wells are described
% using a Peacemann model, giving an extra set of equations that need to be
% assembled. For simplicity, all wells are assumed to be vertical and are
% assigned using the logical (i,j) subindex.

% Set vertical injectors, completed in the lowest 12 layers.
nz = G.cartDims(3);
I = [ 9, 26,  8, 25, 35, 10];
J = [14, 14, 35, 35, 68, 75];
R = [ 4,  4,  4,  4,  4,  4]*1000*meter^3/day;
nIW = 1:numel(I); W = [];
for i = 1 : numel(I)
   W = verticalWell(W, G, rock, I(i), J(i), nz-11:nz, 'Type', 'rate', ...
                    'InnerProduct', 'ip_tpf', ...
                    'Val', R(i), 'Radius', 0.1, 'Comp_i', [1, 0], ...
                    'name', ['I', int2str(i)], 'refDepth', 2500);
end

% Set vertical producers, completed in the upper 14 layers
I = [17, 12, 25, 35, 15];
J = [23, 51, 51, 95, 94];
P = [300, 300, 300, 200, 200];
nPW = (1:numel(I))+max(nIW);
for i = 1 : numel(I)
   W = verticalWell(W, G, rock, I(i), J(i), 1:14, 'Type', 'bhp', ...
                    'InnerProduct', 'ip_tpf', ...
                    'Val', 300*barsa(), 'Radius', 0.1, 'refDepth', 2500, ...
                    'name', ['P', int2str(i)], 'Comp_i', [0, 1]);
end


% Plot grid outline and the wells
figure
myview = struct(...
    'view',     [30,50],      ...  % view angle
    'zoom',     1.8,          ...  % zoom
    'dataasp',  [15 15 2],    ...  % data aspect ratio
    'cb',      'horiz',       ...  % colorbar location
    'dolly',   [0 .3 0],      ...  % camdolly
    'wargs',   {{'height', 30, 'Color', 'k', 'FontSize', 10}},  ...  
    'pargs',   {{'FaceColor', 'none', 'EdgeAlpha', .05}} ...
    );
plotReservoirModel(G, [], W, rmfield(myview,'cb'));
plotGrid(G, vertcat(W(nIW).cells), 'FaceColor', 'b');
plotGrid(G, vertcat(W(nPW).cells), 'FaceColor', 'r');


%% Transmissibilities and initial state
% Initialize solution structures and compute transmissibilities from
% input grid, rock properties, and well structure.
hT = computeTrans(G, rock, 'Verbose', true);
tmult = computeTranMult(G, grdecl);
hT = hT.*tmult;
rSol  = initState(G, W, 0, [0, 1]);

%% Fluid model
fluid      = initSimpleFluid('mu' , [   1,   5]*centi*poise     , ...
                             'rho', [1014, 859]*kilogram/meter^3, ...
                             'n'  , [   2,   2]);


%% Prepare plotting of saturations
figure
myview.zoom = 2.05; myview.dolly = [0,.2,0];
h = plotReservoirModel(G, [], W, myview);
colormap(flipud(winter))
set(h,'Position',[.13 .07 .77 .05]);
[hs,ha] = deal([]); caxis([0 1]);
drawnow

%% Main loop
% In the main loop, we alternate between solving the transport and the flow
% equations. The transport equation is solved using the standard implicit
% single-point upwind scheme with a simple Newton-Raphson nonlinear solver.
T    = 12*year();
dplt = .25*year;
pv   = poreVolume(G,rock);
t    = 0;  plotNo = 1;
dt   = diff(0:year/12:T);
dt   = [dt(1).*sort(repmat(2.^-[1:5 5],1,1)) dt(2:end)];
wSol = cell(numel(dt),1);
for n=1:numel(dt)
   
   rSol = incompTPFA(rSol, G, hT, fluid, 'wells', W);
   rSol = implicitTransport(rSol, G, dt(n), rock, fluid, 'wells', W);

   % Check for inconsistent saturations
   assert(max(rSol.s(:,1)) < 1+eps && min(rSol.s(:,1)) > -eps);

   % Extract well data
   wSol{n} = getWellSol(W, rSol, fluid);
   
   % Increase time and continue if we do not want to plot saturations
   t = t + dt(n);
   if ( t < plotNo*dplt && t <T), continue, end

   % Plot saturation
   delete([hs, ha])
   hs = plotCellData(G, rSol.s(:,1), find(rSol.s(:,1) > 0.01),'EdgeColor','none');
   ha = annotation('textbox', [0 0.93 0.35 0.07], 'String', ...
       sprintf('Water saturation at %5.2f years', convertTo(t,year)));
   drawnow
   plotNo = plotNo+1;
end
