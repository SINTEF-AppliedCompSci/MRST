% testTopSurfaceGrid
clear
gravity on
require ad-fi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Quantities to vary:                        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = 40; %100; %100;
L = 2000; %5000;
dim3 = 1;
H = 15;
dy = 1000;
total_time=5*year; %500*year;
nsteps=3;
dt = total_time/nsteps; % total_time/1000;
perm = 100;
K=100*milli*darcy();
phi= 0.1;

depth=1000
p_press=200;
rate=(H*phi*L*dy)*0.0*day/year;
% set resolution in z-direction
%sigma=dy*h0*sqrt(pi)/(volume/phi);


theta = -0.0*pi/180; %1*pi/180;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Make grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%dx     = L/(n+1); dy = 1000;
%x_grid = [[n/2:-1:1]*(-dx),0,[1:n/2]*dx];
%z_grid = sin(theta).*x_grid;
%z_grid(:) = z_grid(end:-1:1);
clear deck;
cartdims = [n 1 dim3];
nc=prod(cartdims);
% define runspec
deck.RUNSPEC.cartDims=cartdims;
deck.RUNSPEC.DIMENS=cartdims;
deck.RUNSPEC.OIL=1;
%only for eclipse
deck.RUNSPEC.WATER=1;
%deck.RUNSPEC.FIELD=1;
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
chop = @(x) min(max(0,x),1);
s_wc = 0.0; %0.2;
s_or = 0.0; %0.2;
s = linspace(s_wc, 1 - s_or, 2)'; alpha = 2;
s_star = (s - s_wc)/(1 - s_wc - s_or);

deck.PROPS.DENSITY = [600 1000 1];
drho=(deck.PROPS.DENSITY(2)-deck.PROPS.DENSITY(1));
%swof = [s, chop(s_star.^alpha), chop((1-s_star).^alpha), (1-s).0.0];
swof = [s, chop(s_star.^alpha), chop((1-s_star).^alpha), ((1-s)*drho*norm(gravity)*H)/barsa];
%swof = [swof; [1.0 1.0 0.0 0.0]];
deck.PROPS.SWOF{1} = swof;



%deck.PROPS.DENSITY=[900 1000 0.044000000000000];
deck.PROPS.DENSITY = [600 1000 1];
deck.PROPS.ROCK=[100 0 NaN NaN NaN NaN];
%deck.PROPS.ROCK=[4000 3.000000000000000e-06 NaN NaN NaN NaN];
deck.PROPS.PVTW=[100 1.0 0.0 0.40 0];
deck.PROPS.PVCDO=[100 1.0 0.0 0.1 0];
%deck.PROPS.PVDO{1} = [ [300, 800, 8000]'*barsa, [1.05, 1.02, 1.01]', [2.85, 2.99, 3]' ];
%deck.PROPS.PVDO{1}=[100 1 0.1;1000 1 0.1];
% define summary
deck.SUMMARY=[];
% difine SC
%%
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
deck.GRID=grdeclSloping([n, 1, dim3],[L dy H],'theta',theta,'amp',H/5,'lambda',L/2);
deck.GRID.ZCORN=deck.GRID.ZCORN+min(deck.GRID.ZCORN)+depth;
deck.GRID.ACTNUM=int32(ones(nc,1));

deck.GRID.PORO=ones(nc,1)*phi;
deck.GRID.PERMX=ones(nc,1)*perm;
deck.GRID.PERMY=ones(nc,1)*perm;
deck.GRID.PERMZ=ones(nc,1)*perm;

grdecl = grdeclSloping([n, 1, dim3],[L dy H],'theta',theta,'amp',H/5,'lambda',L/4);
g = processGRDECL(grdecl);
g = computeGeometry(g);
g_plot = g; g_plot = computeGeometry(g_plot);
maxx=max(g.nodes.coords(:,1));
g_top = topSurfaceGrid(g);
g_top.cells.H=H*ones(g_top.cells.num,1);
g_top.columns.dz=ones(numel(g_top.columns.cells),1)*H/dim3;
g_top.columns.z = cumulativeHeight(g_top);
%
deck.SOLUTION.SWAT=ones(nc,1);
deck.SOLUTION.SOIL=zeros(nc,1);
% define needed quantitites for simulation
%mu=[deck.PROPS.PVDO{1}(1,3),deck.PROPS.PVTW(1,4)]*centi*poise();
%rho=deck.PROPS.DENSITY(1:2);
deck_c=convertDeckUnits(deck);
G=initEclipseGrid(deck_c)
G=computeGeometry(G);
deck.SOLUTION.PRESSURE=p_press*barsa+(norm(gravity)*deck.PROPS.DENSITY(2))*G.cells.centroids(:,3);

deck_dir='test_deck_well'
writeDeck(deck,deck_dir)

% %% compute timestep
deck = convertDeckUnits(deck);

fluid_old= initEclipseFluid(deck)
G = initEclipseGrid(deck);
G = computeGeometry(G);
state0 = initEclipseState(G, deck, fluid_old);

state0.pressure=p_press*barsa+(G.cells.centroids(:,3)-G.cells.centroids(end,3))*norm(gravity)*deck.PROPS.DENSITY(2);


figure(1),clf
plotGrid(G),view([0 -1 0])

rock  = initEclipseRock(deck);
rock  = compressRock(rock, G.cells.indexMap);

fluid = initEclipseADIFluid(deck);
schedule = deck.SCHEDULE;
systemOW =      initADISystem({'Oil', 'Water'}, G, rock, fluid);

[wellSols      states]      = runScheduleADI(state0, G, rock, systemOW, schedule);
%%

g_top=topSurfaceGrid(G);
%%
figure(2),clf
for nn=1:numel(states)
    clf
    state=states{nn};
    subplot(2,2,1),cla
    title('pressure')
    plotCellData(G,state.pressure);colorbar
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
    plot(g_top.cells.centroids(:,1),state.pressure);
    drawnow;
    pause(0.01)
end
