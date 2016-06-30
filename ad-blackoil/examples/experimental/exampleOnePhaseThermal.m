%% Read the problem from a deckfile

mrstModule add ad-props  ad-core ad-blackoil

% Create grid
%G=cartGrid([100 10 1],[1000 100 10])
G=cartGrid([100 100],[4000 300]);

% Set up the rock structure
rock.perm  = 1000*milli*ones(G.cells.num,1)*darcy;
rock.poro  = ones(G.cells.num,1)*0.1;
rock.lambdaR=ones(G.cells.num,1)*4;
% Create fluid
fluid = initSimpleADIFluid('mu', [1 0.1 1], 'rho', [1 1 1], 'n', [2 2 2]);
fluid.relPerm =@(sW) deal(fluid.krW(sW),fluid.krO(1-sW));
fluid.pvMultR  =@(p) 1+1e-5*(p-200*barsa)/barsa; %to avoid well closing due to incomp.
% Oil rel-perm from 2p OW system.
% Needed by equation implementation function 'eqsfiOWExplictWells'.
fluid.krO = fluid.krOW;

% Get schedule


% Enable this to get convergence reports when solving schedules
verbose = false;

%% Compute constants
% Once we are happy with the grid and rock setup, we compute
% transmissibilities. For this we first need the centroids.
G = computeGeometry(G);
T = computeTrans(G, rock);
p_res=200*barsa;
T_bhp=350;
T_res=300;
%% define wells and timesteps
if(G.griddim==3)
    W = verticalWell([], G, rock,  1,   1, (1:G.cartDims(3)),     ...
        'Type', 'bhp', 'Val', 100*barsa+p_res, ...
        'Radius', 0.4, 'Name', 'P1','Comp_i',[0 1],'sign',1);
    
    W = verticalWell(W, G,rock,  G.cartDims(1),  G.cartDims(2), (1:G.cartDims(3)),     ...
        'Type', 'bhp', 'Val', p_res-1*barsa, ...
        'Radius', 0.4, 'Name', 'I1','Comp_i',[1 0],'sign',-1);
else
    dims=floor(G.cartDims/2);
    wc=sub2ind(G.cartDims,dims(1),1);
    W = addWell([], G, rock,  wc,     ...
        'Type', 'bhp', 'Val', 100*barsa+p_res, ...
        'Radius', 1, 'Name', 'P1','Comp_i',[0 1],'sign',1);
    %{
    W = addWell(W, G,rock,  G.cells.num,     ...
        'Type', 'bhp', 'Val', p_res-1*barsa, ...
        'Radius', 0.4, 'Name', 'I1','Comp_i',[1 0],'sign',-1);
    %}
end
    
dt=diff(linspace(0,400*10,10)*day);
W_c={W};
step=struct('control',ones(numel(dt),1),'val',dt);
schedule=struct('control',struct('W',W_c),'step',step);

% Set up reservoir
% We turn on gravity and set up reservoir and scaling factors.
gravity on

clear wModel
clear nonlinear

clear state;
Cv=4.2e3;
rock.cR=ones;
rock.rhoR=1000;
fluid.bW=@(p,T) 1+(p-p_res)*1e-4/barsa;
fluid.hW=@(p,T) Cv*T;
fluid.uW=@(p,T) Cv*T;
state.pressure = ones(G.cells.num,1)*p_res;
state.s = repmat([1 0],G.cells.num,1);
state.T = ones(G.cells.num,1)*T_res;
%state.wellSols= initWellSolLocal(W, state);
grav=zeros(1,G.griddim);grav(G.griddim)=10;
wModel = WaterThermalModel(G, rock, fluid,'gravity',grav);%, 'deck', deck);
state.wellSols= initWellSolAD(W, wModel, state);
state=rmfield(state,'s');


bc=pside([],G,'Right',p_res,'sat',1);
bc=pside(bc,G,'Left',p_res,'sat',1);
bc.hW=ones(size(bc.face)).*fluid.hW(p_res,T_res);
bcT=addBCT([],vertcat(bc.face),'temperature',T_res);
%bc=pside(bc,G,'Back',p_res);
for i=1:numel(schedule.control)
    schedule.control(i).bc = bc;
    schedule.control(i).bcT=bcT;
    for j=1:numel(schedule.control.W)       
        bhp=schedule.control(i).W(j).val;
        schedule.control(i).W(j).compi=[1];
        schedule.control(i).W(j).hW=fluid.hW(bhp,T_bhp);
    end
end

[wellSols, states] = simulateScheduleAD(state, wModel, schedule);
%%
figure(1),clf
for i=1:numel(states)
    clf,
    plotCellData(G,states{i}.pressure/barsa);colorbar;caxis([200 300])
    %hold on;plot(G.cells.centroids(wc,1),G.cells.centroids(wc,2),'o','Color','r','MarkerSize',10)
    pause(0.1)
    
end
%%
