% testTopSurfaceGrid
clear
close all
gravity on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%% Quantities to vary:                        %%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%n = 250; %100; %100;
%L = 20000; %5000;
%dim3 = 10;
H = 15;
%dy = 1000;

total_time=30000*year; %500*year;
injection_time= total_time; %/4;
dt = total_time/100; % total_time/1000;
% flux-function
n_flux = 1;
K=(100e-3)*darcy();
phi= 0.1;

% set resolution in z-direction
%sigma=dy*h0*sqrt(pi)/(volume/phi);

theta = 0*pi/180; %-1*pi/180;



%% Make grid   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if false
if (true)
   dim3=1;
   g = cartGrid([50, 1, dim3], [50e3 5e3 H]);
   g.nodes.coords(:,3)=g.nodes.coords(:,1)*tan(theta)+g.nodes.coords(:,3)+1000;
   g = computeGeometry(g);
   g3D = g;
elseif false
   dim3=1;nn=30;
   Lx=10e3;Ly=10e3;
   g = cartGrid([nn, nn, dim3], [Lx Ly H]);
   nnn=20;
  % g.nodes.coords(:,3)=g.nodes.coords(:,1)*tan(theta)+g.nodes.coords(:,3)+1000;
   g.nodes.coords(:,3)=sin(nnn*g.nodes.coords(:,1)/Lx).*sin(nnn*g.nodes.coords(:,2)/Ly)*H+g.nodes.coords(:,2)*tan(theta)+g.nodes.coords(:,3)+200;
   g = computeGeometry(g);
end   

grdecl = [];
g_top = topSurfaceGrid(g, 'grdecl', grdecl);
g_top.cells.H=H*ones(g_top.cells.num,1);
g_top.columns.dz=ones(numel(g_top.columns.cells),1)*H/dim3;
g_top.columns.z = cumulativeHeight(g_top);

X=reshape(g_top.cells.centroids(:,1),g_top.cartDims);
Y=reshape(g_top.cells.centroids(:,2),g_top.cartDims);


else
   grdecl = simpleGrdecl([10, 1, 10], 0, 'undisturbed', true);
   
   g3D = processGRDECL(grdecl);
   g3D = computeGeometry(g3D);
    
   [g, g_top, region3D] = make2D3Dgrid([10, 1, 10], [50e3 5e3 H], region3D);
   


end

g_plot = g; g_plot = computeGeometry(g_plot);
maxx=max(g.nodes.coords(:,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Rock and fluid data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rho=[600 1000];mu=[0.4 0.4]*centi*poise;sr=0.0;sw=0.0;
rock.perm  = K*ones(g3D.cells.num,1);
rock.poro  = phi*ones(g3D.cells.num,1);

% rockCoupled.perm = rock.perm(g.cells.indexMap,:);
% rockCoupled.poro = rock.perm(g.cells.indexMap,:);


rock2d.perm  = K*ones(g_top.cells.num,1);
rock2d.poro  = phi*ones(g_top.cells.num,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize fluid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h_init = zeros(size(g_top.cells.H));
%
%h_init(2:9) = H-3;
%zeros(n,1)+min(h0*exp(-sigma*(g_top.cells.centroids(:,1)-L/3).^2),H);
s_init= h_init*(1-sw)/H;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define well and boundary for the test case
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rate = 1e5/year/10;
i = any(g_top.faces.neighbors==0, 2);  % find all outer faces
I = i(g_top.cells.faces(:,1));         % vector of all faces of all cells, true if outer
j = false(4,1);                     % mask, cells can at most have 6 faces,
%j(1:2)=true;  %   extract east, west, north, south
j(1:2)=true;%   extract east
J = j(g_top.cells.faces(:,2));         % vector of faces per cell, true if E,W,N,S
bc_faces = g_top.cells.faces(I & J, 1);
 %depth * density * gravity   
bc_h = addBC([], bc_faces, 'pressure', g_top.faces.z(bc_faces)*1000*norm(gravity));
bc_h = rmfield(bc_h,'sat');
bc_h.h=zeros(size(bc_h.face));             


i = any(g.faces.neighbors==0, 2);  % find all outer faces
I = i(g.cells.faces(:,1));         % vector of all faces of all cells, true if outer
j = false(6,1);                     % mask, cells can at most have 6 faces,
%j(1:2)=true;  %   extract east, west, north, south
j(1:2)=true;%   extract east
J = j(g.cells.faces(:,2));         % vector of faces per cell, true if E,W,N,S
bc_faces = g.cells.faces(I & J, 1);
                                
bc_s = addBC([], bc_faces, 'pressure',g_top.faces.z(bc_faces)*1000*norm(gravity), 'sat', [0 1]);

bc_3D = addBC([], bc_faces, 'pressure', (g3D.faces.centroids(bc_faces,3)+H/2)*1000*norm(gravity), 'sat', [0 1]);



%}
% place well in middle of domain
ind=floor(g_top.cartDims/2);
ind_well=double(g_top.cells.indexMap(sub2ind(g_top.cartDims,ind(1),max(ind(2),1))));
W_h = addWell([], g_top, rock2d, ind_well,'Type', 'rate','Val',rate, ...
               'Radius',0.1,'comp_i',[1 0]);
W_s = [];
W_s = addWell(W_s, g_top, rock2d, ind_well,'Type', 'rate','Val',rate, ...
             'Radius',0.1,'comp_i',[1 0]);
          
          
% W_coupled = addWell([], g, rockCoupl, ind_well,'Type', 'rate','Val',rate, ...
%              'Radius',0.1,'comp_i',[1 0]);          
%           
W_3D = addWell([], g, rock, ind(1),'Type', 'rate','Val',rate,...
               'Radius',0.1,'comp_i',[1 0]);
%W_3D = addWell(W_3D, g, rock, g_top.cells.num/,'Type', 'bhp','Val',100*barsa,'Radius',0.1,'comp_i',[1 0]);

%W_3D = W_s;
%W_s = [];
%W_3D = [];

%W_h = W_s;

src_s = []; %addSource([], ind_well, rate,'sat',[1 0]);
src_h = src_s;

for i=1:numel(src_h)
src_h.sat = nan;
src_h.h = g_top.cells.H(src_s.cell);
end

for i=1:numel(W_h)
   %W_s(i).compi=W_s(i).compi(1:2);
   W_h(i).compi=nan;
   W_h(i).h = g_top.cells.H(W_h(i).cells);
end


%% Define regions for the grid, true for 3D region
region3D = false(g.cells.num,1);
region3D(ind_well-4:ind_well+4) = true;
%region3D(2:9) = true;
g = addGtopFields(g, g_top, g, region3D);

%% Define solvers to solve the system with different discretisations
%

% 1) original version in h formulation
%%{
n=1;
problem{n}.fluid      = initVEFluid(g_top, 'mu', mu, 'rho', rho,'sr',sr,'sw',sw);
problem{n}.sol        = initResSolVE(g_top, 0);
S=computeMimeticIPVE(g_top,rock2d);
%S_2d = computeMimeticIPVE(g_top, rock2D,'Innerproduct','ip_simple');
problem{n}.S  = S;
problem{n}.W=W_h;
problem{n}.bc=bc_h;
problem{n}.src=src_h;
problem{n}.psolver =@(sol,fluid,W,bc, src)...
   solveIncompFlowVE(sol,g_top,S,rock2d,fluid,'wells',W,'bc',bc, 'src', src );
problem{n}.tsolver =@(sol,fluid,dt,W,bc, src)...
   explicitTransportVE(sol, g_top, dt, rock, fluid, ...
   'computeDt', true, 'intVert_poro', false,'dt_dynamic',false,'intVert',false,'wells',W,'bc',bc, 'src', src);
problem{n}.compute_h=false;
problem{n}.compute_sat=true;
problem{n}.col='r';

%  2) Mimetic and implict with s formulation on 3D grid and not integration in vertical
%     direction
%}
n=2;
problem{n}.fluid = initSimpleVEFluid('mu' , mu , 'rho', rho, ...
                           'height'  , g_top.cells.H,...
                           'sr', [sr, sw]);                        
problem{n}.sol        = initResSolVE(g, 0);
S = computeMimeticIP(g, rock); %T=computeTrans(g_top,rock2d);
problem{n}.S=S;
problem{n}.bc=bc_s;
problem{n}.src=src_s;

problem{n}.psolver =@(sol,fluid,W,bc, src)...
   solveIncompFlow(sol, g, S, fluid,'wells',W,'bc',bc, 'src', src);
problem{n}.W=W_3D;
problem{n}.tsolver =@(sol,fluid,dt,W,bc, src)...
         implicitTransport(sol, g, dt, rock, fluid,'vert_avrg',true, ...
                             'wells',W,'bc',bc, 'src', src);
problem{n}.compute_h=true;
problem{n}.compute_sat=false;
problem{n}.col='gs';


n=3;
%  2) TPFA and coupled 3D grid and not integration in vertical
%     direction


problem{n} = problem{n-1};

problem{n}.fluid = initCoupledVEFluid('mu' , mu , 'rho', rho, ...
                           'height'  , g_top.cells.H,...
                           'sr', [sr, sw], 'region3D', region3D, 'n', [1, 1], 'g', g);  
T = computeTrans(g, rock);   
roblem{n}.W=W_3D;
problem{n}.psolver =@(sol,fluid,W,bc, src)...
                     incompTPFAVE_coupled(sol, g, T, fluid,'wells', W, ...
                                    'bc',bc, 'src', src, 'region3D', region3D);
problem{n}.tsolver =@(sol,fluid,dt,W,bc, src)...
implicitTransportVE_coupled(sol, g, dt, rock, fluid,'vert_avrg',true, ...
                            'wells',W,'bc',bc, 'src', src, 'region3D', region3D, 'verbose', true);

problem{n}.col='b*';



n=4;
%  2) TPFA with 3D grid 

% problem{n}.fluid = initCoreyFluid('mu' , mu , 'rho', rho, ...
%                                  'sr', [sr, sw],'kwm', [1, 1], 'n', [1 1]); 
                              
     problem{n}.fluid = initSimpleVEFluid('mu' , mu , 'rho', rho, ...
                           'height'  , g_top.cells.H,...
                           'sr', [sr, sw]);                               
problem{n}.sol        = initResSolVE(g, 0);
T = computeTrans(g, rock);    
problem{n}.psolver =@(sol,fluid,W,bc, src)...
                     incompTPFA(sol, g, T, fluid,'wells', W, ...
                                    'bc',bc, 'src', src);
problem{n}.tsolver =@(sol,fluid,dt,W,bc, src)...
implicitTransport(sol, g, dt, rock, fluid, ...
                            'wells',W,'bc',bc, 'src', src);
                             
problem{n}.bc=bc_3D;
problem{n}.src=src_s;

problem{n}.W=W_3D;
problem{n}.compute_h=true;
problem{n}.compute_sat=false;
problem{n}.col='c--';


%{
n=4;
problem{n} = problem{n-1};                
problem{n}.psolver =@(sol,fluid,W,bc, src)...
   solveIncompFlowVE_coupled(sol, g, S, fluid,'wells',W,'bc',bc, ...
                             'src', src, 'region3D', region3D);
problem{n}.col='c--';



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
                                             


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% initialize all problems with init solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for kk=1:numel(problem)
   problem{kk}.sol.s=s_init;
   problem{kk}.sol.h=h_init;
   problem{kk}.sol.s_max=s_init;
   problem{kk}.sol.max_h=h_init;
   problem{kk}.p_time=0;
   problem{kk}.t_time=0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run transport-simulation:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1);
t = 0;
while t < total_time
   figure(1),clf
   hold on 
   for kk=1:numel(problem)
     disp(kk) 
        
     %if(int_vert_poro)
%          pv_3D(g_top.columns.cells)=rock.poro(g_top.columns.cells)...
%             .*rldecode(g_top.cells.volumes,diff(g_top.cells.columnPos));
%          free_volume = sum(integrateVertically(pv_3D, problem{kk}.sol.h, g_top, cum_z));
      %else
         free_volume= sum(problem{kk}.sol.h.*rock2d.poro.*g_top.cells.volumes);
      %end
      disp(['Total volume free ', num2str(free_volume/1e6)]);
     %{
      if(isfield(sol,'max_h'))
         if(int_vert_poro)
            %trapped_volume= sum(integrateVertically(pv_3D, problem{kk}.sol.max_h-problem{kk}.sol.h, g_top, cum_z))*fluid.sr;
         else
            trapped_volume=sum((problem{kk}.sol.max_h-sol.h).*rock2d.poro.*g_top.cells.volumes)*fluid.sr;
         end
         disp(['Total volume trapped ', num2str(trapped_volume/1e6)]);
         disp(['Total volume ', num2str((trapped_volume+free_volume)/1e6)]);
         times = [times,t];
         volumes = [volumes;[trapped_volume/1e6,free_volume/1e6]];
      end
      %} 
          
      
    
      if(t<injection_time)
         W=problem{kk}.W;
         bc=problem{kk}.bc;
         src = problem{kk}.src;
      %{ 
         if kk == 4, 
               W=problem{1}.W;
         bc=problem{1}.bc;
         end
      %}
      else
         W=[];
         bc=[];
         src = [];
      end
      %% pressure solve
      tmp=tic;
      problem{kk}.sol = problem{kk}.psolver(problem{kk}.sol,problem{kk}.fluid,W,bc, src);
      tmp=toc(tmp);
      problem{kk}.p_time=problem{kk}.p_time+tmp;
      %% transport solve
      tmp=tic;
      nn=1;    
       for ii=1:nn  
         
         if(t>injection_time)
           %problem{kk}.sol.flux(:) = 0;
        
         end
     
      %flux = rand(g.cartDims(1)*2+2,1)*4.8883e-04;        
       %  problem{kk}.sol.flux(:)=0;
        % problem{kk}.sol.flux(1:numel(flux)) = flux;
         problem{kk}.sol = problem{kk}.tsolver(problem{kk}.sol,problem{kk}.fluid,dt/nn,W,bc, src);      
       end      
       
      tmp=toc(tmp);
      problem{kk}.t_time=problem{kk}.t_time+tmp;
      
      if(problem{kk}.compute_h)
         if kk == 3            
             h = sat2height(problem{kk}.sol.s, g_top, rock);
             problem{kk}.sol.h = h;
             s = height2Sat(problem{kk}.sol, g_top, problem{kk}.fluid);
             problem{kk}.sol.s3D = s; 
             
             problem{kk}.sol.s3D(g.cells.inx3D(region3D)) =  problem{kk}.sol.s(region3D); 
         else
         
         
         [h,h_max]=problem{kk}.fluid.sat2height(problem{kk}.sol);
         problem{kk}.sol.h = h;
         problem{kk}.sol.max_h = h_max;
         end
      end
      if(problem{kk}.compute_sat)         
         problem{kk}.sol.s =  (problem{kk}.sol.h*(1-sw)+(problem{kk}.sol.max_h-problem{kk}.sol.h)*sr)./g_top.cells.H;
         problem{kk}.sol.s_max =  problem{kk}.sol.max_h*(1-sw)./g_top.cells.H;
      end      
         
      if(true) %g_top.cartDims(2)==1)
         subplot(2, 2, 1)
         hold on 
         plot(g_top.cells.centroids(:,1), problem{kk}.sol.h,problem{kk}.col)
         if(kk==1);axh=axis();else axis(axh);end
         ylim([0 H])
         title('Height')
         subplot(2, 2, 2)
         hold on 
         plot(g_top.cells.centroids(:,1), problem{kk}.sol.max_h,problem{kk}.col)
         if(kk==1);axh_max=axis();else axis(axh_max);end
         ylim([0 H])
         title('Max Height')
         
         subplot(2, 2, 3)
         hold on 
         plot(g_top.cells.centroids(:,1), problem{kk}.sol.pressure,problem{kk}.col)
         %if(kk==1); axs=axis()else axis(axs);end         
         %axis([min(g_top.cells.centroids(:,1)),max(g_top.cells.centroids(:
         %,1)),0,1])
         title('Pressure')
         subplot(2, 2, 4)
         hold on 
         plot(g_top.cells.centroids(:,1), problem{kk}.sol.flux(1:g.cells.num),problem{kk}.col)
         %if(kk==1); axs_max=axis();else axis(axs_max);end         
         %axis([min(g_top.cells.centroids(:,1)),max(g_top.cells.centroids(:,1)),0,1])
         title('flux')
         %{
         subplot(2, 2, 3)
         hold on 
         plot(g_top.cells.centroids(:,1), problem{kk}.sol.s(:,1),problem{kk}.col)
         if(kk==1);axs=axis();else axis(axs);end         
         ylim([0 1])
         %axis([min(g_top.cells.centroids(:,1)),max(g_top.cells.centroids(:
         %,1)),0,1])
         title('Saturation')
         subplot(2, 2, 4)
         hold on 
         plot(g_top.cells.centroids(:,1), problem{kk}.sol.s_max(:,1),problem{kk}.col)
         if(kk==1);axs_max=axis();else axis(axs_max);end         
         ylim([0 1])
         %axis([min(g_top.cells.centroids(:,1)),max(g_top.cells.centroids(:,1)),0,1])
         title('Max Saturation')
         %}         
      else         
         subplot(numel(problem),2,(2*(kk-1))+1)
         pcolor(X,Y,reshape(problem{kk}.sol.h,g_top.cartDims))
         if(kk==1)cxs=caxis();else caxis(cxs);end;
         title(['Height ',num2str(kk)]);colorbar,shading interp
         subplot(numel(problem),2,(2*(kk-1))+2)
         pcolor(X,Y,reshape(problem{kk}.sol.max_h,g_top.cartDims))
         caxis(cxs)
         title(['Max Height',num2str(kk)]);colorbar,shading interp
      end
         
   end
   t=t+dt;
end


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



