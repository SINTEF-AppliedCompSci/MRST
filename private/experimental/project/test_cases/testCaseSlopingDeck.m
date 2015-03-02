% testTopSurfaceGrid
clear
gravity on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%% Quantities to vary:                        %%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = 20; %100; %100;
L = 2000; %5000;
dim3 = 1;
H = 15;
dy = 1000;
volume = 3e7; h0=30; 
total_time=1*year; %500*year;
nsteps=1;
dt = total_time/nsteps; % total_time/1000;
% flux-function
n_flux = 1;
perm = 100;
K=100*milli*darcy();
phi= 0.1;

% set resolution in z-direction
sigma=dy*h0*sqrt(pi)/(volume/phi);

% refine most near top in z
refine = true; 
theta = -1*pi/180; %1*pi/180;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Make grid   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%dx     = L/(n+1); dy = 1000;
%x_grid = [[n/2:-1:1]*(-dx),0,[1:n/2]*dx];
%z_grid = sin(theta).*x_grid;
%z_grid(:) = z_grid(end:-1:1);
use_deck=true;
if(use_deck)
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
    h_init = zeros(n,1)+min(h0*exp(-sigma*(g_top.cells.centroids(:,1)-L/3).^2),H);
    sr=0.3;sw=0.3;
    s_init= h_init*(1-sw)/H;
    deck.SOLUTION.SWAT=ones(nc,1).*(1-s_init);
    deck.SOLUTION.SOIL=ones(nc,1).*s_init;
    % define needed quantitites for simulation
    mu=[deck.PROPS.PVDO{1}(1,3),deck.PROPS.PVTW(1,4)]*centi*poise();
    rho=deck.PROPS.DENSITY(1:2);
else
    g = cartGrid([n, 1, dim3], [L dy H]);
    g.nodes.coords(:,3)=g.nodes.coords(:,1)*tan(theta)+g.nodes.coords(:,3)+H/5*sin(2*pi*g.nodes.coords(:,1)/(L/4));
    g = computeGeometry(g);
    rho=[600 1000];mu=[0.1 0.4]*centi*poise;sr=0.3;sw=0.3;
    drho=400;
    h_init = zeros(n,1)+min(h0*exp(-sigma*(g_top.cells.centroids(:,1)-L/3).^2),H);
    s_init= h_init*(1-sw)/H;
    g_plot = g; g_plot = computeGeometry(g_plot);
    maxx=max(g.nodes.coords(:,1));
    g_top = topSurfaceGrid(g);
    g_top.cells.H=H*ones(g_top.cells.num,1);
    g_top.columns.dz=ones(numel(g_top.columns.cells),1)*H/dim3;
    g_top.columns.z = cumulativeHeight(g_top);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Rock and fluid data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



rock.perm  = K*ones(g.cells.num,1);
rock.poro  = phi*ones(g.cells.num,1);
rock2d.perm  = K*ones(g_top.cells.num,1);
rock2d.poro  = phi*ones(g_top.cells.num,1);
% the fluid in mrst has been changed: need to multiply by centipoise


%% Initialize solution structure
n=0;
%%{

n=n+1;
problem{n}.fluid      = initVEFluid(g_top, 'mu', mu, 'rho', rho,'sr',sr,'sw',sw);
problem{n}.sol        = initResSol(g_top, 0);
problem{n}.tsolver =@(sol,fluid,dt) explicitTransportVE(sol, g_top, dt, rock, fluid, ...
    'computeDt', true, 'intVert_poro', false,'intVert',false);
problem{n}.compute_h=false;
problem{n}.compute_sat=true;
problem{n}.col='r';
%}
n=n+1;
problem{n}.fluid = initSimpleVEFluid('mu' , mu , 'rho', rho, ...
    'height'  , g_top.cells.H,...
    'sr', [sr, sw]);
problem{n}.sol        = initResSolVE(g_top, 0);
%problem{n}.tsolver =@(sol,fluid,dt) implicitTransportVE(sol, g_top, dt, rock2d, fluid,'vert_avrg',true);
problem{n}.tsolver =@(sol,fluid,dt) implicitTransport(sol, g_top, dt, rock2d, fluid);
problem{n}.compute_h=true;
problem{n}.compute_sat=false;
problem{n}.col='b*';
%{
n=n+1;
problem{n}=problem{n-1};
problem{n}.fluid = initVEVerticalIntegratedFluid('mu' , mu, ...
                                                'rho', rho, ...
                                                'G'  , g_top,...
                                                'sr', [sr, sw],...
                                                'rockORG',rock);
problem{n}.sol        = initResSolVE(g_top, 0);
problem{n}.compute_h=true;
problem{n}.compute_sat=false;
%problem{n}.tsolver =@(sol,fluid,dt) implicitTransportVE(sol, g_top, dt, rock2d, fluid,'vert_avrg',true);
problem{n}.tsolver =@(sol,fluid,dt) implicitTransport(sol, g_top, dt, rock2d, fluid);
problem{n}.col='gx';
%}
% generate all solvers for h formulation
%{
intvertPoro_vec=[true,false]
dt_dynamic_vec=[true,false]
intVert_vec=[true,false]
for intvertPoro=intvertPoro_vec
   for dt_dynamic=dt_dynamic_vec
      for intVert=intVert_vec
         n=n+1
         problem{n}=problem{1}
         problem{n}.tsolver =@(sol,fluid,dt) explicitTransportVE(sol, g_top, dt, rock, fluid, ...
            'computeDt', true, 'intVert_poro', intvertPoro,'dt_dynamic',dt_dynamic,'intVert',intVert);
      end
   end
end
%}




%%
for kk=1:numel(problem)
    problem{kk}.sol.s=s_init;
    problem{kk}.sol.h=h_init;
    problem{kk}.sol.extSat=[s_init,s_init];
    problem{kk}.sol.h_max=h_init;
    %problem{kk}.sol=rmfield(problem{kk}.sol,'max_h');
    %problem{kk}.sol=rmfield(problem{kk}.sol,'s_max');
end
%% Run transport-simulation:
t = 0;
h1=figure(1)
h2=figure(2)
while t < total_time
    %figure(1),clf
    %figure(2),clf
    hold on
    for kk=1:numel(problem)
        problem{kk}.sol = problem{kk}.tsolver(problem{kk}.sol,problem{kk}.fluid,dt);
        if(problem{kk}.compute_h)
            [h,h_max]=problem{kk}.fluid.sat2height(problem{kk}.sol);
            problem{kk}.sol.h = h;
            problem{kk}.sol.h_max = h_max;
            problem{kk}.sol.s_max = problem{kk}.sol.extSat(:,2);
        end
        if(problem{kk}.compute_sat)
            problem{kk}.sol.s =  (problem{kk}.sol.h*(1-sw)+(problem{kk}.sol.h_max-problem{kk}.sol.h)*sr)./g_top.cells.H;
            problem{kk}.sol.s_max =  problem{kk}.sol.h_max*(1-sw)./g_top.cells.H;
        end
    end
    
    
    set(0, 'CurrentFigure', h1);
    clf
    for kk=1:numel(problem)
        subplot(2, 2, 1)
        hold on
        plot(g_top.cells.centroids(:,1), problem{kk}.sol.h,problem{kk}.col)
        ylim([0 h0*1.3]);
        title('Height')
        subplot(2, 2, 2)
        hold on
        plot(g_top.cells.centroids(:,1), problem{kk}.sol.h_max,problem{kk}.col)
        ylim([0 h0*1.3]);
        title('Max Height')
        subplot(2, 2, 3)
        hold on
        plot(g_top.cells.centroids(:,1), problem{kk}.sol.s(:,1),problem{kk}.col)
        if(t==0)
            ax=axis();
        end
        %axis([min(g_top.cells.centroids(:,1)),max(g_top.cells.centroids(:
        %,1)),0,1])
        title('Saturation')
        subplot(2, 2, 4)
        hold on
        plot(g_top.cells.centroids(:,1), problem{kk}.sol.s_max,problem{kk}.col)
        axis(ax);
        %axis([min(g_top.cells.centroids(:,1)),max(g_top.cells.centroids(:,1)),0,1])
        title('Max Saturation')
    end
    
    set(0, 'CurrentFigure', h2);
    clf
    %%{
    hold on
    plot(g_top.cells.centroids(:,1),g_top.cells.z,'LineWidth',2)
    plot(g_top.cells.centroids(:,1),g_top.cells.z+g_top.cells.H,'LineWidth',2)
    for kk=1:numel(problem)
         plot(g_top.cells.centroids(:,1),g_top.cells.z+problem{kk}.sol.h,problem{kk}.col)
         plot(g_top.cells.centroids(:,1),g_top.cells.z+problem{kk}.sol.h_max,problem{kk}.col)        
    end
    ylim([min(g_top.cells.z),max(g_top.cells.z+g_top.cells.H)]);
    set(gca,'YDir','reverse');
    %}
    drawnow;


t=t+dt;
end
%% write deck and VE simulation result
deck_dir='sloping_case_deck_nowell'
writeDeck(deck,deck_dir)
save(fullfile(deck_dir,'s_init.txt'),'-ascii','s_init')
tmp=problem{1}.sol.s;
save(fullfile(deck_dir,'s_end.txt'),'-ascii','tmp')
% %% compute timestep
% sinth= (z_grid(2:end)-z_grid(1:end-1))'/dx; sinth= [sinth(1);sinth;sinth(end)];
% const_perp = norm(gravity())*K*drho*sqrt(1-sinth.^2);
% const_par = norm(gravity())*K*drho*sinth;
% t_inject=20*year; t = 0; hh=0:H/100:H;
% v=[-well_rate*ones(n_well,1);well_rate*ones(n_well,1)]/(2*dy*H); c2 = 0;
% lam1 = @(h,h_res) h/fluid.mu(1); lam2 = @(h,h_res) ((H-h)+h_res*c2)/fluid.mu(2);
% FF = @(h,h_res) lam1(h,h_res).*lam2(h,h_res)./(lam1(h,h_res)+lam2(h,h_res));
% dd_max=max(abs( (FF(hh(1:end-1),0)-FF(hh(2:end),0))./(hh(1:end-1)-hh(2:end))));
% dt_perp=min(0.1./(2*pi*abs(const_perp)*dd_max*max(H)./(dx^2)));
% dt_par = min(0.5./(abs(const_par)*dd_max/(dx*phi)));
% dt_adv=min(0.1./(well_rate/(dx*dy*H)));
% %NB:
% %dt = abs(min([dt_perp,dt_par,dt_adv]))
%



