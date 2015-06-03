clearvars -except METISPATH mrstVerbose screenSize
min_size = 1; cell_size = 5;
a = 1/25; % fracture aperture
dof_frac = 1; % Coarse dof
dof_matrix = 50; % Coarse dof
layers = 5;
%% User defined fracture lines and appropriate grid

nx = 100; ny = 100;
G = cartGrid([nx ny]);
G = computeGeometry(G);
fl = [30 30 70 70; 30 70 70 30];
% load('Lines13_Network1_100by100Domain.mat'); nx = G.cartDims(1); ny = G.cartDims(2);
flayers = 2:4;


%% Process fracture lines into grid

[G,fracture] = processFracture2D(G,fl); fracture.aperture = a;
% figure; plotFractureLines(G,fracture,'network'); % 3rd arg can be line #'s, network or none
% figure; plotMarkedCells(G, fracture,'lines'); % 3rd arg can be line #'s, network or none


%% Compute CI and grid fracture

dispif(mrstVerbose, 'Computing CI and gridding fracture...\n\n');
G = CIcalculator2D(G,fracture);
% plotCI(G,fracture);

[G,F,fracture] = gridFracture2D(G,fracture,'min_size',min_size,'cell_size',cell_size);
% figure; plotFractureNodes2D(G,F,fracture,show_node_num); clear show_node_num


%% Make Layered Grid
dispif(mrstVerbose, 'Extruding...\n\n');
Gl = makeLayers(G,flayers,layers);


%% Rock and fluid properties
Gl.rock.perm = ones(Gl.cells.num,1)*darcy();
Gl.rock.poro = 0.2*ones(Gl.cells.num, 1);
K_frac = 10000; % Scaling factor = K_frac/K_mat D
poro_frac = 0.5;
Gl = makeRockFrac(Gl, K_frac, 'permtype','homogeneous','rockporo',poro_frac);
fluid = initSimpleFluid('mu' , [   1,  1] .* centi*poise     , ...
    'rho', [1000, 700] .* kilogram/meter^3, ...
    'n'  , [   1,   1]);
s = linspace(0,1,101)'; kr=fluid.relperm(s);


%% Global grid with NNC's and corresponding transmissibility
dispif(mrstVerbose, 'Define NNC''s...\n\n');
[Gl,T] = makeNNCgrid(G,Gl,F,fracture,flayers);
G = Gl; clear Gl


%% Init

pinit = 5*barsa();
state = initResSol (G, pinit);

radius = 1e-2;
cellinj = 1:nx*ny:nx*ny*layers;
cellprod = nx*ny:nx*ny:nx*ny*layers;
W   = addWell([], G.Matrix, G.Matrix.rock, cellinj, 'Type', 'rate',...
    'Val', 0.1/day(), 'Radius', radius, 'Sign',1, 'Comp_i', [1, 0]);
W   = addWell(W, G.Matrix, G.Matrix.rock, cellprod, ...
    'Type', 'bhp', 'Val', 0.5*barsa(), 'Radius', radius, 'Sign', -1, 'Comp_i', [0, 1]);

bc = [];
state = incompTPFA(state, G, T, fluid,  ...
    'bc',bc, 'Wells', W, 'MatrixOutput', true, 'use_trans',true);


%% Saturation
dispif(mrstVerbose, 'Transport...\n\n');
pv     = poreVolume(G,G.rock);
nt     = 30;
t90   = 0.9*(sum(pv)/sum(state.wellSol(1).flux));
Time   = t90;
dT     = Time/nt;
dTplot = Time/3;  % plot only every 1500 days
N      = fix(Time/dTplot);
t = 0;
pvi = zeros(nt,1); s = zeros(G.cells.num,nt);
count = 0;
screenSize = get(0,'screensize');
figure('Position',[screenSize(1:end-1),screenSize(end)/1.5]);
plotNo = 1;
while t < Time,
    count = count+1;
    state = implicitTransport(state, G, dT, G.rock, fluid, 'wells', W, 'bc', bc, 'Trans', T);

    % Check for inconsistent saturations
    s(:,count) = state.s;
    assert(max(s(:,count)) < 1+eps && min(s(:,count)) > -eps);

    % Update solution of pressure equation.
    state  = incompTPFA(state, G, T, fluid, 'wells', W, 'bc', bc, 'use_trans',true);
    
    t = t + dT;
    if ( t < plotNo*dTplot && t <Time), continue, end
    % Plot saturation
    pvi(count) = 100*(sum(state.wellSol(1).flux)*t)/sum(pv);
    heading = [num2str(round(pvi(count),1)),  ' % PVI'];
    subplot(1,N,plotNo);
    plotWell(G,W); view(-10,40);
    plotCellData(G,state.s,unique(G.nnc.cells(:,2)));
    plotCellData(G,state.s,setdiff(1:G.cells.num,unique(G.nnc.cells(:,2))),'FaceAlpha',0.5,'EdgeAlpha',0);
%     plotCellData(G,state.s,'FaceAlpha',0.7,'EdgeAlpha',0);
    colormap(flipud(gray))
    title(heading);
    
    plotNo = plotNo+1;
end
figure;
plotToolbar(G,s,'FaceAlpha',0.7,'EdgeAlpha',0.01);
colormap(flipud(gray))
