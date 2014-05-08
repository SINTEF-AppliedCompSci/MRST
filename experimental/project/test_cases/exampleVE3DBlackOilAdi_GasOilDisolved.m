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
close all
mrstVerbose true
gravity on
[nx,ny,nz] = deal(50, 1, 1);     % Cells in Cartsian grid
[Lx,Ly,H]  = deal(2000,1000,50); % Physical dimensions of reservoir
total_time = 5*year;             % Total simulation time
nsteps     = 40;                 % Number of time steps in simulation
dt         = total_time/nsteps;  % Time step length
perm       = 100;                % Permeability in milli darcies
phi        = 0.1;                % Porosity
depth      = 1000;               % Initial depth
ipress     = 300;
dp         = 70;

%% Create input deck and construct grid
% Create an input deck that can be used together with the fully-implicit
% solver from the 'ad-fi' module. Since the grid is constructed as part of
% setting up the input deck, we obtain it directly. 
[deck, G] = sinusDeckAdi_GasOilDisolved([nx ny nz], [Lx Ly H], nsteps, dt, ...
                         -.1*pi/180, depth, phi, perm, ...
                         (H*phi*Lx*Ly)*0.2*day/year, ipress,dp);

% Alternatively, we could read deck from file and construct the grid
% deck = readEclipseDeck( ...
%    fullfile(VEROOTDIR,'data','decks','sinusDeckAdi.DATA');
% G = initEclipseGrid(deck);

figure, plotGrid(G),view([0 -1 0]), box on
                      


%% Initialize data structures
% First, we convert the input deck to SI units, which is the unit system
% used by MRST. Second, we initialize the rock parameters from the deck;
% the resulting data structure may have to be post-processed to remove
% inactive cells. Then we set up the fluid object and tell the ad-fi solver
% that that we are working with an oil-gas system.
deck  = convertDeckUnits(deck);
rock  = initEclipseRock(deck);
rock  = compressRock(rock, G.cells.indexMap);
fluid = initDeckADIFluid(deck);
% set the capillary pressure and the VE relperms explicitely
Gt = topSurfaceGrid(G);
explicit_VE=true
relperm_caase='hystersis';

if(explicit_VE)
    pressure_case='simple_time_mix';
    %pressure_case='simple_instant';
    switch pressure_case
        case 'simple_time_mix'
            l_fac=1e-10
            fluid.bO=@(po,rs,flag,varargin) (po-200*barsa)*l_fac+1;
            fluid.BO=@(po,rs,flag,varargin) 1./fluid.bO(po,rs,flag,varargin);
            fluid.dis_rate=5e-13;
        case 'simple_instant'
            l_fac=1e-6
            fluid.bO=@(po,rs,flag,varargin) (po-200*barsa)*l_fac+1;
            fluid.BO=@(po,rs,flag,varargin) 1./fluid.bO(po,rs,flag,varargin);
        otherwise
    end
    
    
    fluid.bG=@(pg,varargin) pg*0.0+1;
    fluid.BG=@(pg,varargin) pg*0.0+1;
    fluid.muO=@(po,rs,flag,varargin) 0.4e-3*(po*0+1);
    fluid.muG=@(pg,varargin) 1e-4*(pg*0+1);
    dis_max=0.01;
    fluid.rsSat=@(po,rs,flag,varargin)   (po*0+1)*dis_max;
    
    relperm_case='hysteresis';
    %relperm_case='simple';
    switch relperm_case
        case 'simple'
            fluid.krG=@(sg,varargin) sg;
            fluid.krOG=@(so,varargin) so;
            fluid.pcOG=@(sg,varargin) norm(gravity)*(fluid.rhoOS-fluid.rhoGS)*(sg).*Gt.cells.H;
            fluid=rmfield(fluid,'relPerm');
            %fluid.re
        case 'hysteresis'
            res_gas = 0.2;
            fluid = addVE3DRelperm(fluid,'res_oil',0,'res_gas',res_gas,'Gt',Gt);            
        otherwise
            disp('Use deck as fluid')                
    end
end

systemOG  = initADISystem({'Oil', 'Gas','DisGas'}, G, rock, fluid)%, 'tol_mb', 0);
%systemOW.nonlinear.linesearch=true;
if(strcmp(pressure_case,'simple_instant'))
    systemOG.nonlinear.maxIterations=10;
else
    systemOG.nonlinear.maxIterations=20;
end
%% Run the schedule setup in the file
% Before we can run the schedule, we make sure that we have an initial
% hydrostatic pressure distribution. Then we pick the schedule from the
% input deck and start the simulator.
%x0 = initEclipseState(G, deck, initEclipseFluid(deck));
z  = G.cells.centroids(:,3);
clear x0;
x0.pressure = ipress*barsa +(z(:)-z(end))*norm(gravity)*deck.PROPS.DENSITY(2);
x0.s(:,1)=deck.SOLUTION.SOIL;
x0.s(:,2)=deck.SOLUTION.SGAS;
x0.rs=ones(size(deck.SOLUTION.RS))*dis_max*0.0;
x0.smax=x0.s;
x0.smin=x0.s;
%x0.s=x0.s(:,[2,1]);
%x0.s=x0.z(:,[2,1]);


[wellSols, states] = runScheduleADI(x0, G, rock, systemOG, deck.SCHEDULE);

%% Plot results
%figure

xc = Gt.cells.centroids(:,1);
zt = Gt.cells.z;
zb = zt + Gt.cells.H;
for nn=1:numel(states)
    clf
    state=states{nn};
    rsH=state.s(:,1).*state.rs/dis_max;
    subplot(2,2,1),cla
    title('pressure')
    plotCellData(G,state.pressure/barsa);colorbar('horiz'), caxis([100 600])
    
    subplot(2,2,2),cla
    title('saturation')
    plotCellData(G,state.s(:,2));colorbar('horiz'), caxis([0 1])
    
    % plot as VE
    subplot(2,2,3),cla
   
   
    plot(xc,state.pressure/barsa); set(gca,'YLim',[100 600]);

    subplot(2,2,4),cla,hold on
    %%
     %figure()
     clf
    patch(xc([1 1:end end]), [zt(end)-20; zt; zt(end)-20],.7*[1 1 1]);
    patch(xc([1 1:end end]), [zb(end)+20; zb; zb(end)+20],.7*[1 1 1]);
    patch(xc([1:end end:-1:1]), ...
      [zt + Gt.cells.H.*state.s(:,2); zt(end:-1:1)], 'g')
    patch(xc([1:end end:-1:1]), ...
      [zt + Gt.cells.H.*state.s(:,2); zt(end:-1:1)+Gt.cells.H.*(state.smax(end:-1:1,2))],[1 1 0])
    patch(xc([1:end end:-1:1]), ...
      [zt + Gt.cells.H.*state.smax(:,2); zt(end:-1:1)+Gt.cells.H.*(state.s(end:-1:1,2)+rsH(end:-1:1))], [0 0.1 1])
    patch(xc([1:end end:-1:1]), ...
      [zt + Gt.cells.H.*(state.s(:,2)+rsH); zb(end:-1:1)], [0 0 0.6])
  
  
    line(xc([1:end]),zt + Gt.cells.H.*state.smax(:,2),'LineWidth',2,'Color','k')
    set(gca,'YDir','reverse'), axis tight
    %%
    drawnow;
    pause(0.1)
end
