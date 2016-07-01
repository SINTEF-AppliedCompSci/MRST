%% Example not using deckformat and only using simple one phase solver in 
%  the object oriented AD framework
%  The example only in its basic 2D version only have one well and boundary
%  condtions.

mrstModule add ad-props  ad-core ad-blackoil

% Create grid
G=cartGrid([100 100],[1000 100]);

% Set up the rock structure
rock.perm  = 100*milli*ones(G.cells.num,1)*darcy;
rock.poro  = ones(G.cells.num,1)*0.3;

% Create fluid
%fluid = initSimpleADIFluid('mu', [1 0.1 1], 'rho', [1 1 1], 'n', [2 2 2]);
fluid = initSimpleADIFluid('mu', 1*centi*poise, 'rho', 1000,'phases','W');
fluid.pvMultR  =@(p) 1+1e-5*(p-p_res)/barsa; %to avoid well closing due to incomp.
fluid.bW=@(p) 1+(p-p_res)*1e-4/barsa;

% Enable this to get convergence reports when solving schedules
verbose = false;

%% Compute constants
% Once we are happy with the grid and rock setup, we compute
% transmissibilities. For this we first need the centroids.
G = computeGeometry(G);
T = computeTrans(G, rock);
p_res=200*barsa;
%% define wells and timesteps
if(G.griddim==3)
    W = verticalWell([], G, rock,  1,   1, (1:G.cartDims(3)),     ...
        'Type', 'bhp', 'Val', 100*barsa+p_res, ...
        'Radius', 0.4, 'Name', 'P1','Comp_i',1,'sign',1);
else
    dims=floor(G.cartDims/2);
    wc=sub2ind(G.cartDims,dims(1),1);
    W = addWell([], G, rock,  wc,     ...
        'Type', 'bhp', 'Val', 100*barsa+p_res, ...
        'Radius', 0.1, 'Name', 'P1','Comp_i',1,'sign',1);
end


% Set up reservoir
% We turn on gravity and set up reservoir and scaling factors.
gravity on
clear wModel
clear nonlinear
clear state;

state.pressure = ones(G.cells.num,1)*p_res;
grav=zeros(1,G.griddim);grav(G.griddim)=-10;

% Initiate model
wModel = WaterModel(G, rock, fluid,'gravity',grav);

%% Equilibrate model with gravity without well and boundary condtions
% define initail bc for equilibration
if(G.griddim==2)
    bc=pside([],G,'Front',p_res,'sat',1);
else
    bc=pside([],G,'Top',p_res,'sat',1);
end
state.wellSol=initWellSolAD([], wModel, state);
nonlinear=NonLinearSolver;
[state, status] = nonlinear.solveTimestep(state, 10000*day, wModel,'bc',bc)
clf,plotCellData(G,state.pressure/barsa),colorbar


%%
% Define boundary, conditions, and wells
% set time steps
dt=diff(linspace(0,1,20)*day);
W_c={W};
step=struct('control',ones(numel(dt),1),'val',dt);
schedule=struct('control',struct('W',W_c),'step',step);
for i=1:numel(schedule.control)
    schedule.control(i).bc = [];%bc;
    for j=1:numel(schedule.control.W)
        schedule.control(i).W(j).compi=[1];
    end
end

% Run simulation
[wellSols, states] = simulateScheduleAD(state, wModel, schedule);

%% Plot simulation
%
figure(1),clf
for i=1:numel(states)
    clf,
    plotCellData(G,states{i}.pressure/barsa);colorbar;%caxis([200 300])
    hold on;plot(G.cells.centroids(wc,1),G.cells.centroids(wc,2),'o','Color','r','MarkerSize',10)
    pause(0.1)
    
end
%%
