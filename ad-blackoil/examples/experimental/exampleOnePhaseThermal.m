%% Read the problem from a deckfile

mrstModule add ad-props  ad-core ad-blackoil

% Create grid
% G=cartGrid([100 10 1],[1000 100 10])
G=cartGrid([100 100],[1000 100]);
% For faster 2D plotting
G.cells.sortedCellNodes=getSortedCellNodes(G);

% Set up the rock structure
rock.perm  = 100*milli*ones(G.cells.num,1)*darcy;
rock.poro  = ones(G.cells.num,1)*0.1;
rock.lambdaR=ones(G.cells.num,1)*4;
% Create fluid
p_res=200*barsa;
fluid = initSimpleADIFluid('mu', 1*centi*poise, 'rho', 1000,'phases','W');
fluid.pvMultR  =@(p) 1+1e-5*(p-p_res)/barsa; %to avoid well closing due to incomp.
fluid.bW=@(p) 1+(p-p_res)*1e-4/barsa;
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
        'Radius', 0.3, 'Name', 'P1','Comp_i',[0 1],'sign',1);
end
for i=1:numel(W)
	W(i).lims=inf;
end
dt=diff(linspace(0,100,20)*day);
W_c={W};
step=struct('control',ones(numel(dt),1),'val',dt);
schedule=struct('control',struct('W',W_c),'step',step);

% Set up reservoir
% We turn on gravity and set up reservoir and scaling factors.
gravity on

clear wModel
clear nonlinear
clear state;

%% define simple thermal properties and bW  and possibly dependent on termperature
%NB: uW and hW not thermodynamically consistent.
Cv=4.2e3;
rock.cR=ones;
rock.rhoR=1000;
fluid.bW=@(p,T) 1+(p-p_res)*1e-4/barsa;
fluid.muW=@(p,T) fluid.muW(p);
fluid.hW=@(p,T) Cv*T;
fluid.uW=@(p,T) Cv*T;


% define initial state
state.pressure = ones(G.cells.num,1)*p_res;
state.T = ones(G.cells.num,1)*T_res;

% define gravity zero gravity
grav=zeros(1,G.griddim);%grav(G.griddim)=-10;

% make model
wModel = WaterThermalModel(G, rock, fluid,'gravity',grav);%, 'deck', deck);


% define boundary conditions
bc = pside([],G,'Right',p_res,'sat',1);
bc = pside(bc,G,'Left',p_res,'sat',1);
bc.hW = ones(size(bc.face)).*fluid.hW(p_res,T_res);
bcT = addBCT([],vertcat(bc.face),'temperature',T_res);
for i=1:numel(schedule.control)
    schedule.control(i).bc = bc;
    schedule.control(i).bcT=bcT;
    for j=1:numel(schedule.control.W)       
        bhp=schedule.control(i).W(j).val;
        schedule.control(i).W(j).compi=[1];
        schedule.control(i).W(j).hW=fluid.hW(bhp,T_bhp);
    end
end

% solve  system
[wellSols, states] = simulateScheduleAD(state, wModel, schedule);
%%
figure(1),clf
for i=1:numel(states)
    clf,
    subplot(2,1,1),cla
    plotCellData(G,states{i}.pressure/barsa);colorbar;caxis([200 300]) % plot of pressure
    subplot(2,1,2),cla
    plotCellData(G,states{i}.T);colorbar;caxis([300 350]) % plot of temperature
    %hold on;plot(G.cells.centroids(wc,1),G.cells.centroids(wc,2),'o','Color','r','MarkerSize',10)
    pause(0.1)
    
end

%% Copyright notice

% <html>
% <p><font size="-1">
% Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.
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
