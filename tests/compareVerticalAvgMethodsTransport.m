% testTopSurfaceGrid
clear
gravity on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%% Quantities to vary:                        %%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = 20; %100; %100;
L = 2000; %5000;
dim3 = 10;
H = 15;
dy = 1000;
volume = 3e7; h0=30; 
total_time=40*year; %500*year;
dt = total_time/10; % total_time/1000;
% flux-function
n_flux = 1;
K=(100e-3)*darcy();
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
    grdecl = grdeclSloping([n, 1, dim3],[L dy H],'theta',theta,'amp',H/5,'lambda',L/4);
    g = processGRDECL(grdecl);
    g = computeGeometry(g);
else
    g = cartGrid([n, 1, dim3], [L dy H]);
    g.nodes.coords(:,3)=g.nodes.coords(:,1)*tan(theta)+g.nodes.coords(:,3)+H/5*sin(2*pi*g.nodes.coords(:,1)/(L/4));
    g = computeGeometry(g);
end

g_plot = g; g_plot = computeGeometry(g_plot);
maxx=max(g.nodes.coords(:,1));
g_top = topSurfaceGrid(g);
g_top.cells.H=H*ones(g_top.cells.num,1);
g_top.columns.dz=ones(numel(g_top.columns.cells),1)*H/dim3;
g_top.columns.z = cumulativeHeight(g_top);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Rock and fluid data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


rho=[600 1000];mu=[0.1 0.4]*centi*poise;sr=0.3;sw=0.3;
rock.perm  = K*ones(g.cells.num,1);
rock.poro  = phi*ones(g.cells.num,1);
rock2d.perm  = K*ones(g_top.cells.num,1);
rock2d.poro  = phi*ones(g_top.cells.num,1);
% the fluid in mrst has been changed: need to multiply by centipoise
drho=400;
h_init = zeros(n,1)+min(h0*exp(-sigma*(g_top.cells.centroids(:,1)-L/3).^2),H);
s_init= h_init*(1-sw)/H;

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






