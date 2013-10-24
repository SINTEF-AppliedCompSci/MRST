%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Example to test vertical average caculation using different formulsation
% The test case uses two wells at the start
% We compare to simulations
%   1) formulation using h as variable and mimetic for pressure and
%      explicit time stepping. This was the original formulation in the
%      vertical average module.
%   2) Using s as variable and using tpfa an implicit transport from the
%      mrst core. 
%
%   the example also set up a compate deck file which can be written and
%   may be used for simulation by tradition solvers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
gravity on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%% Quantities to vary:                        %%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = 30; %100; %100;
L = 2000; %5000;
dim3 = 1; % cells in the vertical direction
H = 15; % hight of resevoir
dy = 1000; % with of slize
total_time=100*year; %500*year;
nsteps=20; % number of steps
dt = total_time/nsteps; % total_time/1000;
injection_time=total_time/10; % injection time
%n_flux = 1; 
perm = 100;
K=perm*milli*darcy();
phi= 0.1;

% initial related quantities
depth=1000;  % depth of reservoir
p_press=200; % initial pressure
rate=(H*phi*L*dy)*0.2*day/injection_time;


% refine most near top in z
theta = -1*pi/180; %1*pi/180;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Make a deck struct   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    %deck.PROPS.SWOF{1}=[s,s.^alpha,(1-s).^alpha,s*0];
    chop = @(x) min(max(0,x),1);
    s_wc = 0.0; %0.2;
    s_or = 0.0; %0.2;
    s = linspace(s_wc, 1 - s_or, 2)'; alpha = 2;
    s_star = (s - s_wc)/(1 - s_wc - s_or);
    swof = [s, chop(s_star.^alpha), chop((1-s_star).^alpha), s*0.0];
    %swof = [s, chop(s_star.^alpha), chop((1-s_star).^alpha), s*drho*norm(g)*0.0];
    %swof = [swof; [1.0 1.0 0.0 0.0]];
    deck.PROPS.SWOF{1} = swof;
    
    pres = convertTo(convertFrom(6000, psia), barsa);
    deck.SOLUTION.PRESSURE=ones(nc,1)*pres;
     
    %deck.PROPS.DENSITY=[900 1000 0.044000000000000];
    deck.PROPS.DENSITY = [600 1000 1];
    deck.PROPS.ROCK=[100 1.0e-6 NaN NaN NaN NaN];
    %deck.PROPS.ROCK=[4000 3.000000000000000e-06 NaN NaN NaN NaN];
    deck.PROPS.PVTW=[100 1.0 0.0 0.40 0];
    %deck.PROPS.PVDO{1} = [ [300, 800, 8000]'*barsa, [1.05, 1.02, 1.01]', [2.85, 2.99, 3]' ];
    deck.PROPS.PVDO{1}=[100 1 0.1;1000 1 0.1];
    % define summary
    deck.SUMMARY=[];
    % difine SC
    %%
    % define grid
    deck.SCHEDULE.step.val=ones(nsteps,1)*dt/day;
    %deck.SCHEDULE.step.control=ones(nsteps,1);
    deck.GRID=grdeclSloping([n, 1, dim3],[L dy H],'theta',theta,'amp',H/5,'lambda',L/4);
    deck.GRID.ZCORN=deck.GRID.ZCORN+depth;
    deck.GRID.ACTNUM=int32(ones(nc,1));
    
    deck.GRID.PORO=ones(nc,1)*phi;
    deck.GRID.PERMX=ones(nc,1)*perm;
    deck.GRID.PERMY=ones(nc,1)*perm;
    deck.GRID.PERMZ=ones(nc,1)*perm;
    
    grdecl = grdeclSloping([n, 1, dim3],[L dy H],'theta',theta,'amp',H/5,'lambda',L/4);
    
    %
   
    sr=0.3;sw=0.3;
    deck.SOLUTION.SWAT=ones(nc,1);
    deck.SOLUTION.SOIL=ones(nc,1).*0.0;
    % define needed quantitites for simulation
    mu=[deck.PROPS.PVDO{1}(1,3),deck.PROPS.PVTW(1,4)]*centi*poise();
    rho=deck.PROPS.DENSITY(1:2);
    deck.SCHEDULE.control.WELSPECS=...
    {...
    'I01'    'W'    [  ceil(cartdims(1)/2)]    [ ceil(cartdims(2)/2) ]    [1000]    'OIL'    [0]    'STD'    'SHUT'    'NO'    [0]    'SEG'    [0];...
    'P01'    'W'    [  cartdims(1)]    [cartdims(2)]                     [1004.03647]    'OIL'    [0]    'STD'    'SHUT'    'NO'    [0]    'SEG'    [0];...
    };%#ok
    radius=0.01;
    deck.SCHEDULE.control.COMPDAT=...
    {...
    'I01'     [  ceil(cartdims(1)/4)]    [ ceil(cartdims(2)/2) ]   [1]    [cartdims(3)]    'OPEN'    [0]    [0]    [radius]    [-1]    [0]    'Default'    'Z'    [-1];...
    'P01'    [  cartdims(1)]    [cartdims(2)]     [1]    [cartdims(3)]                                         'OPEN'    [0]    [0]    [radius]    [-1]    [0]    'Default'    'Z'    [-1];...
    };%#ok

% use scaled spe10 rates
deck.SCHEDULE.control.WCONINJE=...
    {...
    'I01'  'OIL'  'OPEN'  'RESV'  [rate]  [rate]  [300-89.4018]  [Inf]  [0]  [0]...
    };%#ok
deck.SCHEDULE.control.WCONPROD=...
    {...
    'P01'  'OPEN'  'BHP'  [Inf]  [Inf]  [Inf]  [Inf]  [Inf]  [p_press-89.4018]  [0]  [0]  [0];...
    };%#ok
deck.SCHEDULE.control=[deck.SCHEDULE.control;deck.SCHEDULE.control];
deck.SCHEDULE.control(2).WCONINJE{3}='SHUT';
% define grid
deck.SCHEDULE.step.val=ones(nsteps,1)*dt/day;
deck.SCHEDULE.step.val=[deck.SCHEDULE.step.val,deck.SCHEDULE.step.val*40];

%% deck is finnished
%writeDeck('test_slope',deck);
% convert units
deck_converted=convertDeckUnits(deck);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Rock and grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make grid
g = processGRDECL(deck_converted.GRID);% process grid
% and alternative would be g = initEclipseGrid(deck_conveted);
g = computeGeometry(g);% calculate geometry
maxx=max(g.nodes.coords(:,1));
% make top surface grid
g_top = topSurfaceGrid(g);
g_top.cells.H=H*ones(g_top.cells.num,1);
g_top.columns.dz=ones(numel(g_top.columns.cells),1)*H/dim3;
g_top.columns.z = cumulativeHeight(g_top);
% get permeability
rock.perm  = deck_converted.GRID.PERMX;
rock.poro  = deck_converted.GRID.PORO;
rock2d  = averageRock(rock, g_top);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define well and boundary for the test case
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bondary
bc_s=[];
bc_h = bc_s;
%bc_h = rmfield(bc_h,'sat');
%bc_s.sat=zeros(size(bc_h.face));
%bc_h.h=zeros(size(bc_h.face));


%% make well from deck which have to be one cell layers
require deckformat
% process the well information we only use control(1)
W_3D = processWells(g,rock,deck_converted.SCHEDULE.control(1));
W_3D_2ph=W_3D;
% change resv to rate which is used in the 2ph simulations
for i=1:numel(W_3D)    
   if(strcmp(W_3D(i).type,'resv'))
    W_3D_2ph(i).type='rate';
   end
   %  well indices is correct since they are calculated in 3D which 
   % include the hight in the WI
end
% make VE wels for s and h formulation
% we will use mimetic for the h formulation with ip_simple innerproduct
W_h = convertwellsVE(W_3D_2ph, g, g_top, rock2d,'ip_simple');
% we will use tpfa method for the s formulation 
W_s = convertwellsVE_s(W_3D_2ph, g, g_top, rock2d,'ip_tpf');
% correct the definition of input value according to solver
for i=1:numel(W_h)
   W_s(i).compi=[1 0];
   W_h(i).compi=nan;
   W_h(i).h = g_top.cells.H(W_h(i).cells);
end

%% no sorce in this example
src_s = []; %addSource([], ind_well, rate,'sat',[1 0]);
src_h = src_s;
for i=1:numel(src_h)
    src_h.sat = nan;
    src_h.h = g_top.cells.H(src_s.cell);
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define solvers to solve the system with
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% H based solver
n=1;
problem{n}.fluid      = initVEFluidHForm(g_top, 'mu', mu, 'rho', rho,'sr',sr,'sw',sw); %fluid
problem{n}.sol        = initResSol(g_top, 0);% initial state
S=computeMimeticIPVE(g_top,rock2d,'Innerproduct','ip_simple'); % inner product 
problem{n}.W=W_h;% well
problem{n}.bc=bc_h; % boundary condition
problem{n}.src=src_h; % sources
% define pressure solver
problem{n}.psolver =@(sol,fluid,W,bc, src)...
   solveIncompFlowVE(sol,g_top,S,rock2d,fluid,'wells',W,'bc',bc, 'src', src );
% define transport solver
problem{n}.tsolver =@(sol,fluid,dt,W,bc, src)...
   explicitTransportVE(sol, g_top, dt, rock, fluid, ...
   'computeDt', true, 'intVert_poro', false,'intVert',false,'wells',W,'bc',bc, 'src', src);
problem{n}.compute_h=false; % this solver compute h 
problem{n}.compute_sat=true; % this solver do not compute saturation so we have to do that separately to compeare
problem{n}.col='r'; % color for the plotting
% generic twophase solver

n=2;
problem{n}.fluid = initSimpleVEFluid_s('mu' , mu , 'rho', rho, ...
                           'height'  , g_top.cells.H,...
                           'sr', [sr, sw]);                        
problem{n}.sol        = initResSolVE(g_top, 0);
T=computeTrans(g_top,rock2d);
cellno = gridCellNo(g_top);
T=T.*g_top.cells.H(cellno); % modfy the 2D transmissbilities
problem{n}.bc=bc_s;
problem{n}.src=src_s;
problem{n}.psolver =@(sol,fluid, W, bc, src)...
    incompTPFA(sol,g_top,T,fluid,'wells',W,'bc',bc,'src',src);
problem{n}.W=W_s;
problem{n}.tsolver =@(sol,fluid,dt,W,bc, src)...
   implicitTransport(sol, g_top, dt,rock2d, fluid,'wells',W,'bc',bc,'src',src);
problem{n}.compute_h=true;
problem{n}.compute_sat=false;
problem{n}.col='gs';

   


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% initialize all problems with init solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for kk=1:numel(problem)
   s_init=deck.SOLUTION.SOIL;
   h_init=deck.SOLUTION.SOIL.*g_top.cells.H;
   problem{kk}.sol.s=s_init;
   problem{kk}.sol.h=h_init;
   problem{kk}.sol.extSat=[s_init,s_init];
   problem{kk}.sol.h_max=h_init;
   problem{kk}.p_time=0;
   problem{kk}.t_time=0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run transport-simulation:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig1=figure(1);
fig2=figure(2);
fig3=figure(3);
t = 0;
while t < total_time   
   set(0, 'CurrentFigure', fig1);clf
   set(0, 'CurrentFigure', fig2);clf
   set(0, 'CurrentFigure', fig3);clf
   set(gca,'YDir','reverse')
   for kk=1:numel(problem)

      free_volume= sum(problem{kk}.sol.s.*rock2d.poro.*g_top.cells.volumes.*g_top.cells.H);
      disp(['Total volume free ',num2str(kk),' ', num2str(free_volume/1e6)]);

      if(t<injection_time)
         W=problem{kk}.W;
         bc=problem{kk}.bc;
         src = problem{kk}.src;
      else
         W=[];
         bc=[];
         src = [];
      end
      %% pressure solve
      tmp=tic;
      problem{kk}.sol = problem{kk}.psolver(problem{kk}.sol, problem{kk}.fluid,  W, bc, src);
      tmp=toc(tmp);
      problem{kk}.p_time=problem{kk}.p_time+tmp;
      %% transport solve
      tmp=tic;
      nn=1;    
       for ii=1:nn  
         
         if(t>injection_time)
           %problem{kk}.sol.flux(:) = 0;% uncomment if forsing zero total
           %flow
         end
         problem{kk}.sol = problem{kk}.tsolver(problem{kk}.sol,problem{kk}.fluid,dt/nn,W,bc, src);      
      end
      tmp=toc(tmp);
      problem{kk}.t_time=problem{kk}.t_time+tmp;
      
      if(problem{kk}.compute_h)
         % if the solver do not compute h, h_max calculate it
         [h,h_max]=problem{kk}.fluid.sat2height(problem{kk}.sol);
         problem{kk}.sol.h = h;
         problem{kk}.sol.h_max = h_max;
         % set the value s_max for convenience
         problem{kk}.sol.s_max = problem{kk}.sol.extSat(:,2);
      end
      if(problem{kk}.compute_sat)
          % if the solver do not compute s, s_max calculate it
          problem{kk}.sol.s =  (problem{kk}.sol.h*(1-sw)+(problem{kk}.sol.h_max-problem{kk}.sol.h)*sr)./g_top.cells.H;
          problem{kk}.sol.s_max =  problem{kk}.sol.h_max*(1-sw)./g_top.cells.H; 
      end
      
   
      
      
      if(g_top.cartDims(2)==1)
         % plot saturation and max_h
         set(0, 'CurrentFigure', fig1);
         
         subplot(2, 2, 1)
         hold on 
         plot(g_top.cells.centroids(:,1), problem{kk}.sol.h,problem{kk}.col)
         if(kk==1);axh=axis();else axis(axh);end
         ylim([0 H])
         title('Height');xlabel('x'); ylabel('h');
         
         subplot(2, 2, 2)
         hold on
         plot(g_top.cells.centroids(:,1), problem{kk}.sol.h_max,problem{kk}.col)
         if(kk==1);axh_max=axis();else axis(axh_max);end
         ylim([0 H])
         title('Max Height');xlabel('x'); ylabel('h_max') ;  
         
         subplot(2, 2, 3)
         hold on 
         plot(g_top.cells.centroids(:,1), problem{kk}.sol.s(:,1),problem{kk}.col)
         if(kk==1);axs=axis();else axis(axs);end         
         ylim([0 1])
         title('Saturation');xlabel('x'); ylabel('s') 
         
         subplot(2, 2, 4)
         hold on 
         plot(g_top.cells.centroids(:,1), problem{kk}.sol.s_max(:,1),problem{kk}.col)
         if(kk==1);axs_max=axis();else axis(axs_max);end         
         ylim([0 1])
         title('Max Saturation');xlabel('x'); ylabel('s_max')         

         
         
         set(0, 'CurrentFigure', fig2);
         hold on
         plot(g_top.cells.centroids(:,1),problem{kk}.sol.pressure/barsa,problem{kk}.col)
         if(~problem{kk}.compute_sat)
            % also plot co2 pressure the default is to use the persure of
            % the second phase pressure for the incompressible solvers
            plot(g_top.cells.centroids(:,1),(problem{kk}.sol.pressure-problem{kk}.fluid.pc(problem{kk}.sol))/barsa,['-',problem{kk}.col])
         else
            plot(g_top.cells.centroids(:,1),(problem{kk}.sol.pressure-norm(gravity)*deck.PROPS.DENSITY(1)*problem{kk}.sol.h)/barsa,['s',problem{kk}.col]) 
         end
         xlabel('x')
         ylabel('pressure')
         title('Comparing pressure for the different solvers')
         
        %   1) for the injection period the well indexs which is not exact
        %      for a 1D calculation effect the result.
         
         set(0, 'CurrentFigure', fig3);
         hold on
         if(kk==1)
            plot(g_top.cells.centroids(:,1),g_top.cells.z,'k','LineWidth',2)
            plot(g_top.cells.centroids(:,1),g_top.cells.z+g_top.cells.H,'k','LineWidth',2)
            mind=floor(g_top.cells.num/2);
            text(g_top.cells.centroids(mind,1),g_top.cells.z(mind)-5,'Top surface')
            text(g_top.cells.centroids(mind,1),g_top.cells.z(mind)+g_top.cells.H(mind)+5,'Bottom surface')
            mind=floor(g_top.cells.num/4);
            text(g_top.cells.centroids(mind,1)-50,g_top.cells.z(mind)+problem{kk}.sol.h(mind)-2,'Free CO2','Color','b')
            text(g_top.cells.centroids(mind,1),g_top.cells.z(mind)+problem{kk}.sol.h_max(mind)+1,'Max CO2','Color','r') 
         end
         plot(g_top.cells.centroids(:,1),g_top.cells.z+problem{kk}.sol.h,'b')
         plot(g_top.cells.centroids(:,1),g_top.cells.z+problem{kk}.sol.h_max,'r')
         set(gca,'FontSize',16)
         box on;axis tight
         
         title('Surfaces')
         xlabel('x')
         ylabel('depth')
         
         
      else         
         subplot(numel(problem),2,(2*(kk-1))+1)
         pcolor(X,Y,reshape(problem{kk}.sol.h,g_top.cartDims))
         if(kk==1)
             cxs=caxis();
         else
             caxis(cxs);
         end;
         title(['Height ',num2str(kk)]);colorbar,shading interp
         subplot(numel(problem),2,(2*(kk-1))+2)
         pcolor(X,Y,reshape(problem{kk}.sol.max_h,g_top.cartDims))
         caxis(cxs)
         title(['Max Height',num2str(kk)]);colorbar,shading interp
      end
         
   end
   set(0, 'CurrentFigure', fig1);
   drawnow;
   set(0, 'CurrentFigure', fig2);
   drawnow;
   set(0, 'CurrentFigure', fig3);
   
   drawnow;
   t=t+dt;
end
%% the state is now almost stationary we add the analytic calculation of the pressure
set(0, 'CurrentFigure', fig2);
hold on
%%
hold on
%In this example it is always only water at the bottom
%for hydrostatic conditions we can find the pressure at the
%bottom assuming zero pressure at the top of the first column

plot(g_top.cells.centroids(:,1),(1/barsa)*norm(gravity)*deck.PROPS.DENSITY(2)*...
    (g_top.cells.z+g_top.cells.H(1)-(g_top.cells.z(1)))+problem{kk}.sol.pressure(1)/barsa,'dk-')
% water pressure extrapolated hydrostatic to top of reservoir
plot(g_top.cells.centroids(:,1),(1/barsa)*norm(gravity)*deck.PROPS.DENSITY(2)*...
    (g_top.cells.z-(g_top.cells.z(1)))+problem{kk}.sol.pressure(1)/barsa,'dk-')
% water pressure at interface of co2 and water
plot(g_top.cells.centroids(:,1),(1/barsa)*norm(gravity)*deck.PROPS.DENSITY(2)*...
    (g_top.cells.z+problem{1}.sol.h-(g_top.cells.z(1)))+problem{kk}.sol.pressure(1)/barsa,'dk-')

%co2 at the top surface
plot(g_top.cells.centroids(:,1),...
    (1/barsa)*norm(gravity)*deck.PROPS.DENSITY(2)*(g_top.cells.z+g_top.cells.H(1)...
    -(g_top.cells.z(1)))+problem{kk}.sol.pressure(1)/barsa-...
    (1/barsa)*norm(gravity)*(deck.PROPS.DENSITY(2)*(g_top.cells.H-problem{kk}.sol.h)...
    +deck.PROPS.DENSITY(1)*problem{2}.sol.h),'dr-')
%%
text(600,1,'Water pressure at bottom')
text(600,-3.5,'Water extrapolated to top')
text(1300,-2,'Pressure at interface')
text(1850,-3,'Pressure at top')
% In this plot of the pressure at the end we see
%    - The mimetic calculate the pressure at the interface between
%        co2 and water
%      - The tpfa with the given fluid calculate the extrapolated
%        water pressure at the top

