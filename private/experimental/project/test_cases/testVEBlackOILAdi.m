%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Example seting up simple model using standard black oil model to 
%% calculate VE model
%  In this example we show how to put a up a standard format black oil
%  model to calculate VE model. For the actual calculation we use the 
%  solver in MRST based on automatic diffrenciation AD-FI.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
%% assure that gravity is on 
gravity on
% require the ad-fi module
require ad-fi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Quantities to vary:                        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = 40; % number of cells
L = 2000; % length of reservoir
dim3 = 1; % only one cell in the z direction
H = 15; % hight of the reservoir
dy = 1000; % with of the reservoir
total_time=5*year; % total simulation time
nsteps=10; % number of steps
dt = total_time/nsteps; % time steps
perm = 100; %permeability
phi= 0.1; % porosity
% depth of reservoir
depth=1000
p_press=200;
rate=(H*phi*L*dy)*0.2*day/year;

theta = -0.1*pi/180; %1*pi/180;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Make a proper deck struct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear deck;
cartdims = [n 1 dim3];
nc=prod(cartdims);
% define runspec
deck.RUNSPEC.cartDims=cartdims;
deck.RUNSPEC.DIMENS=cartdims;
deck.RUNSPEC.OIL=1;
deck.RUNSPEC.WATER=1;
deck.RUNSPEC.METRIC=1;
deck.RUNSPEC.TABDIMS=[1     1    20    50    20    50     1    20    20     1    10     1    -1     0     1];
deck.RUNSPEC.WELLDIMS=[5 10 2 1 5 10 5 4 3 0 1 1];
deck.RUNSPEC.AQUDIMS=[0 0 0 0 10 10 0 0];
deck.RUNSPEC.START=734139;
% one sat num region
deck.REGIONS.SATNUM=ones(nc,1);
%define props
s=linspace(0,1,2)';alpha=2;
%deck.PROPS.SWOF{1}=[s,s.^alpha,(1-s).^alpha,s*0];
s = linspace(0, 1, 2)'; 
deck.PROPS.DENSITY = [600 1000 1];
drho=(deck.PROPS.DENSITY(2)-deck.PROPS.DENSITY(1));
%% Set up swof with the capillary pressure equal to what it should be in a
% VE model and linear relperm functions
% this work if the top and bottum is parallel if not it is more convenient
% to make black oil fluids which embed the grid and in the fluid object.
% For a standard simulator small modification is also need in the gravity
% term for these cases.
swof = [s, s, 1-s, ((1-s)*drho*norm(gravity)*H)/barsa];
deck.PROPS.SWOF{1} = swof;
deck.PROPS.DENSITY = [600 1000 1];
deck.PROPS.ROCK=[100 0 NaN NaN NaN NaN];
% uncomment if rock compressibility should be included
%deck.PROPS.ROCK=[100 3.000000000000000e-06 NaN NaN NaN NaN];
deck.PROPS.PVTW=[100 1.0 0.0 0.40 0];
deck.PROPS.PVCDO=[100 1.0 0.0 0.1 0];
deck.SUMMARY=[];
%% set up wells
% 
deck.SCHEDULE.control.WELSPECS=...
    {...
    'I01'    'W'    [  ceil(cartdims(1)/2)]    [ ceil(cartdims(2)/2) ]    [1000]    'OIL'    [0]    'STD'    'SHUT'    'NO'    [0]    'SEG'    [0];...
    'P01'    'W'    [  cartdims(1)]    [cartdims(2)]                     [1004.03647]    'OIL'    [0]    'STD'    'SHUT'    'NO'    [0]    'SEG'    [0];...
    };
radius=0.01;
deck.SCHEDULE.control.COMPDAT=...
    {...
    'I01'     [  ceil(cartdims(1)/4)]    [ ceil(cartdims(2)/2) ]   [1]    [cartdims(3)]    'OPEN'    [0]    [0]    [radius]    [-1]    [0]    'Default'    'Z'    [-1];...
    'P01'    [  cartdims(1)]    [cartdims(2)]     [1]    [cartdims(3)]                                         'OPEN'    [0]    [0]    [radius]    [-1]    [0]    'Default'    'Z'    [-1];...
    };
use_pressure=false;

% use scaled spe10 rates
deck.SCHEDULE.control.WCONINJE=...
    {...
    'I01'  'OIL'  'OPEN'  'BHP'  [rate]  [rate]  [300-89.4018]  [Inf]  [0]  [0]...
    };
deck.SCHEDULE.control.WCONPROD=...
    {...
    'P01'  'OPEN'  'BHP'  [Inf]  [Inf]  [Inf]  [Inf]  [Inf]  [p_press-89.4018]  [0]  [0]  [0];...
    };
deck.SCHEDULE.control=[deck.SCHEDULE.control;deck.SCHEDULE.control];
deck.SCHEDULE.control(2).WCONINJE{3}='SHUT';
% define grid
deck.SCHEDULE.step.val=ones(nsteps,1)*dt/day;
deck.SCHEDULE.step.val=[deck.SCHEDULE.step.val,deck.SCHEDULE.step.val*40];
deck.SCHEDULE.step.control=[ones(nsteps,1),2*ones(nsteps,1)];

%% make grid section
deck.GRID=grdeclSloping([n, 1, dim3],[L dy H],'theta',theta,'amp',H/5,'lambda',L/2);
deck.GRID.ZCORN=deck.GRID.ZCORN+min(deck.GRID.ZCORN)+depth;
deck.GRID.ACTNUM=int32(ones(nc,1));

deck.GRID.PORO=ones(nc,1)*phi;
deck.GRID.PERMX=ones(nc,1)*perm;
deck.GRID.PERMY=ones(nc,1)*perm;
deck.GRID.PERMZ=ones(nc,1)*perm;

grdecl = grdeclSloping([n, 1, dim3],[L dy H],'theta',theta,'amp',H/5,'lambda',L/4);

deck.SOLUTION.SWAT=ones(nc,1);
deck.SOLUTION.SOIL=zeros(nc,1);
deck_c=convertDeckUnits(deck);
% we make the grid here to to the initialization of pressure
G=initEclipseGrid(deck_c)
G=computeGeometry(G);
deck.SOLUTION.PRESSURE=p_press/barsa+(norm(gravity)*deck.PROPS.DENSITY(2))*G.cells.centroids(:,3)/barsa;
clear G;
%% write deck to file
deck_dir='test_deck_well'
writeDeck(deck,deck_dir)
%% Start setting up simulation using the deck
%% convert units of deck
% alternatively deck=readEclipseDeck('test_deck_well/test_deck_well.DATA')
deck = convertDeckUnits(deck);
% generate grid 
G = initEclipseGrid(deck);
G = computeGeometry(G);
% make MRST fluid only to use the init original init functions
fluid_old= initEclipseFluid(deck)
state0 = initEclipseState(G, deck, fluid_old);
% Set hydrostatic initial condition explicitly
state0.pressure=p_press*barsa+(G.cells.centroids(:,3)-G.cells.centroids(end,3))*norm(gravity)*deck.PROPS.DENSITY(2);
% alternative state0.pressure=deck.SOLUTION.PRESSURE;

fig1=figure(1),clf
% plot grid
plotGrid(G),view([0 -1 0])
box on
% initialize perm and poro
rock  = initEclipseRock(deck);
% if actum cells rock need to be compressed
rock  = compressRock(rock, G.cells.indexMap);
% init ADI fluid
fluid = initEclipseADIFluid(deck);
% pick schedual section
schedule = deck.SCHEDULE;
% make oil water ADI system
systemOW =      initADISystem({'Oil', 'Water'}, G, rock, fluid);
% run schedual
[wellSols      states]      = runScheduleADI(state0, G, rock, systemOW, schedule);
%%
% make top surface grid to facilitate plotting
g_top=topSurfaceGrid(G);
%%
fig2=figure(2),clf
for nn=1:numel(states)
    clf
    state=states{nn};
    subplot(2,2,1),cla
    title('pressure')
    plotCellData(G,state.pressure/barsa);colorbar
    subplot(2,2,2),cla
    title('saturation')
    plotCellData(G,state.s(:,1));colorbar
    % plot as as VE
    subplot(2,2,3),cla,hold on
    s_co2=state.s(:,2);
    plot(g_top.cells.centroids(:,1),g_top.cells.z,'LineWidth',2)
    plot(g_top.cells.centroids(:,1),g_top.cells.z+g_top.cells.H,'LineWidth',2)
    plot(g_top.cells.centroids(:,1),g_top.cells.z+g_top.cells.H.*s_co2,'r')
    set(gca,'YDir','reverse')
    subplot(2,2,4),cla
    plot(g_top.cells.centroids(:,1),state.pressure/barsa);
    drawnow;
    pause(0.01)
end
