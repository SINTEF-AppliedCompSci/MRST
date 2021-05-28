%% VE simulation in a standard black-oil solver
%  In this example we show how to set up a standard format black-oil
%  model that can be used to simulate a VE model. For the actual
%  simulation,  we use the fully-implicit solver in MRST from the 'ad-fi'
%  module, which is based on automatic differentiation.

try
   require deckformat ad-fi
catch %#ok<CTCH>
   mrstModule add deckformat ad-fi
end

%% Parameters for the simulation
gravity off
[nx,ny,nz] = deal(100, 1, 1);     % Cells in Cartsian grid
[Lx,Ly,H]  = deal(400,10,15); % Physical dimensions of reservoir
total_time = year/10;             % Total simulation time
nsteps     = 50;                 % Number of time steps in simulation
dt         = total_time/nsteps;  % Time step length
perm       = 200;                % Permeability in milli darcies
phi        = 0.3;                % Porosity
depth      = 1000;               % Initial depth
ipress     = 200;                % Initial pressure

% Read and process file.
mycase='simple'
switch mycase
    case 'SPE1'
        fn    = fullfile(ROOTDIR,'\modules\ad-fi\examples\SPE1\', 'odeh_adi.data');
        deck = readEclipseDeck(fn);
        deck = convertDeckUnits(deck);
        fluid = initDeckADIFluid(deck);
    case 'simple'
        fluid=initSimpleADIFluid('mu', [0.4 4 0.04]*centi*poise, 'rho', [1000 700 100], 'n', [2 2 2]);
        fluid.rsSat =@(x) 0*x; 
        fluid.relPerm = @(sw, sg, varargin) deal(fluid.krW(sw, varargin{:}),fluid.krO(1 - sw - sg, varargin{:}),fluid.krG(sg, varargin{:}));


    otherwise
        error();
end

% The deck is given in field units, MRST uses metric.

% Create a special ADI fluid which can produce differentiated fluid
% properties.




% The deck is given in field units, MRST uses metric.
%G = initEclipseGrid(deck);
%G = computeGeometry(G);
%rock  = compressRock(rock, G.cells.indexMap);
G = cartGrid([nx,ny,nz],[Lx,Ly,H])
G = computeGeometry(G);
rock.perm=perm*ones(G.cells.num,3)*milli*darcy;
rock.poro=phi*ones(G.cells.num,1);

% The case includes gravity
gravity off




% for match 1./muO has to be interpolated linearly
%system.getEquations @=  eqsfiBlackOilExplicitWells;
%Wext=processWellsLocal(G, rock, schedule.control(1));
%% add temprature for well related quantities
% add temprature in fluid
cW       = 4.1813*10.^6/1e3; % energy density per mass
cR       = 2.17*10.^6;% energy density per volume
fluid.uW = @(p,T) cW.*T;
fluid.uG = @(p,T) 0.1*cW.*T;
fluid.uO = @(p,T) 0.7*cW.*T;
fluid.uR = @(T) cR.*T;

fluid.hW = @(p,T) cW.*T;
fluid.hG = @(p,T) 0.1*cW.*T;
fluid.hO = @(p,T) 0.7*cW.*T;
system = initADISystem({'Oil','Water','Gas','disgas','T'}, G, rock, fluid, 'cpr', false);
system.nonlinear.tolMB=1e-5
system.nonlinear.tolCNV=1e-7
system.nonlinear.tol=1e-5;
system.nonlinear.use_ecltol=false;

W=[];
W = verticalWell(W, G, rock,  1,   1, (1:G.cartDims(3)),     ...
                     'Type', 'bhp', 'Val', 300*barsa, ...
                     'Radius', 0.125, 'Name', 'I1','Comp_i',[1 0 0]);
W(end).sign=1;
W = verticalWell(W, G, rock,  G.cartDims(1),   G.cartDims(2), (1:G.cartDims(3)),     ...
                     'Type', 'bhp', 'Val', 200*barsa, ...
                     'Radius', 0.125, 'Name', 'P1','Comp_i',[1 0 0]);
W(end).sign=-1;

for i=1:numel(W)
    % black oil relatetd
   W(i).bhpLimit = inf;
   if(W(i).sign>0)
       W(i).lims=struct('bhp',inf','rate',inf,'orat',inf,'wrat',inf,'grat',inf,'lrat',inf)
   else
       W(i).lims=struct('bhp',-inf','rate',-inf,'orat',-inf,'wrat',-inf,'grat',-inf,'lrat',-inf)
   end
   W(i).refDepth = G.cells.centroids(W(i).cells(1),3);
   W(i).dZ=0;
   % termprature related
    W(i).T=400;
    W(i).hW=fluid.hW(200*barsa, W(i).T);
    W(i).hO=fluid.hW(200*barsa, W(i).T);
    W(i).hG=fluid.hW(200*barsa, W(i).T);
end

%calculate rock conductivity
fake_rock.perm=4.0e0*ones(G.cells.num,1);
T = computeTrans(G,fake_rock);
Trans=1./accumarray(G.cells.faces(:,1),1./T,[G.faces.num,1]);
internal=all(G.faces.neighbors>0,2);

system.s.T_r=Trans(internal);

mrst_schedule = struct('step',struct('control',ones(nsteps,1),'val',dt*ones(nsteps,1)));
mrst_schedule.W={W};
system.fluid = fluid;


x0.s=nan(G.cells.num,3);
x0.pressure = ipress*barsa*ones(G.cells.num,1);
x0.s(:,1)=0.05*ones(G.cells.num,1);
x0.s(:,2)=0.8*ones(G.cells.num,1);
x0.s(:,3)=1-sum(x0.s(:,1:2),2);
x0.smax=x0.s;
x0.smin=x0.s;
x0.T=300*ones(G.cells.num,1);
x0.rs = fluid.rsSat(x0.pressure);
x0.rs(x0.s(:,3)==0)=0.4*x0.rs(x0.s(:,3)==0);
x0.rv=0;
state=x0;
%state.T=373*ones(G.cells.num,1);
timer = tic;
myfys='oil'
myfys='boil_temp'
switch myfys
    case 'boil_temp'
        system.stepFunction =@(state0, state, meta, dt, W, G, system, varargin) stepBlackOilTemp(state0, state, meta, dt, G, W, system, fluid);
        system.getEquations =@ eqsfiBlackOilTemp;
        system.updateState  =@  updateStateBlackOilTemp;
        fluid.muW =@(p,T) fluid.muW(p)./(1+1e-2.*(T-300));
        fluid.muO =@(p,rs,isSat,T) fluid.muO(p,rs,isSat)./(1+10e-1.*(T-300));
        system.fluid=fluid;
    case 'oil'

    otherwise
        error()
end

%[wellSols states iter] = runScheduleADI(state, G, rock, system, schedule);
[wellSols states iter schedule] = runMrstADI(state, G, system, mrst_schedule,'force_step',false,'dt_min',0.1*day);
toc(timer)

%% Plot the solution
% We opt for a simple volume plot of the gas saturation. If opengl
% capability is set to software, we fall back to a simpler cell data plot.
% If you have problems with getting good plots you can set useVolume to
% false.
nn=4;
figure(),
subplot(nn,1,1);
xc=G.cells.centroids(:,1);
for i=1:numel(states)
    state=states{i};
    subplot(nn,1,1),cla
    plot(xc,state.pressure/barsa)
    subplot(nn,1,2),cla
    plot(xc,state.s);legend('water','oil','gas')
    subplot(nn,1,3),cla
    plot(xc,state.rs)
    subplot(nn,1,4),cla
    plot(xc,state.T)
    pause(0.01);
end

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
