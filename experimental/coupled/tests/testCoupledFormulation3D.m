% testTopSurfaceGrid
clear
clc
close all
gravity on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%% Quantities to vary:                        %%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% place well in top or bottom of domain
ind=[1 1]; %ceil(g_top.cartDims/2);

% EXAMPLE 1:
%gravity reset on
k_ix = 1; 

% EXAMPLE 2:
gravity(0.1); % reset on
%k_ix = g.cartDims(3);

bc_side = 2;
rate = 1e7/year;
H = 15;

total_time= 500*year;
injection_time= total_time/4;
dt =  total_time/100;
% flux-function for 3D
n_flux = 1;
K=(100e-3)*darcy();
phi= 0.1;

% dip of reservoir
theta = 0*pi/180; %-1*pi/180;

printFigs =false;

%% Make grid   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% grdecl = simpleGrdeclMine([20, 1, 10], [50e3 5e3 H], 0, 'undisturbed', true);
% TUKLET MED
grdecl = simpleGrdecl([20, 5, 10], 0, 'undisturbed', true, 'physDims', [50e3 5e3 H]);
% change simpleGRDECL to take in physdims

%g.nodes.coords(:,3)=g.nodes.coords(:,1)*tan(theta)+g.nodes.coords(:,3)+1000;
g3D = processGRDECL(grdecl);
g3D = computeGeometry(g3D);
clear ijk
[ijk{1:3}] = ind2sub(g3D.cartDims, g3D.cells.indexMap(:));
ijk        = [ijk{:}];
region3D_org = ismember(ijk(:,1), 1:6); %round(g3D.cartDims(1)/2)-4:round(g3D.cartDims(1)/2)+4 );
%region3D = ismember(ijk(:,1), round(g3D.cartDims(1)/2)-4:round(g3D.cartDims(1)/2)+4 ); %& ...
%    ismember(ijk(:,2), round(g3D.cartDims(2)/2)-1:round(g3D.cartDims(2)/2)+1 )  ;
%   region3D = false(g3D.cells.num,1);
%region3D = region3D | (ijk(:,1) == 4 & ijk(:,2) == 1  );

[g, g_top, region3D] = make2D3Dgrid(grdecl, region3D_org);
g_top_coupled = topSurfaceGrid(g);

g_plot = g3D; g_plot = computeGeometry(g_plot);
maxx=max(g.nodes.coords(:,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Rock and fluid data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rho=[600 1000];mu=[0.4 0.4]*centi*poise;sr=0.0;sw=0.0;
rock.perm  = K*ones(g3D.cells.num,1);
rock.poro  = phi*ones(g3D.cells.num,1);

%rock.perm(region3D_org) = rand(nnz(region3D_org),1).*rock.perm(region3D_org);


rockCoupled.perm = rock.perm(g.cells.indexMap,:);
rockCoupled.poro = rock.poro(g.cells.indexMap,:);

rock2d = averageRock(rock, g_top);


rockCoupled.perm(~region3D) = rock2d.perm(g.cells.mapTopSurface(~region3D));

figure;
plotGrid(g, 'faceColor', 'none');
%plotGrid(g, find(region3D));
view([0 0])
axis tight off

if printFigs
set(gcf, 'render', 'painters')
 % shading flat; 
   print -depsc2 fig-g_coupled.eps
end

figure;
plotGrid(g3D, 'faceColor', 'none');
%plotGrid(g, find(region3D));
view([0 0])
axis tight off
if printFigs
set(gcf, 'render', 'painters')
   print -depsc2 fig-g3D.eps
end

figure
%plotGrid(g, 'faceColor', 'none');
plotCellData(g, log10(rockCoupled.perm));
view(3)

clf
plotGrid(g, 'faceColor', 'none');
plotGrid(g, find(region3D));
axis tight off;
view(3)

if printFigs
   print fig-g_region3D.png
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize fluid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h_init = zeros(size(g_top.cells.H));
%
%h_init(2:9) = H-3;
s_init= h_init*(1-sw)/H;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define well and boundary for the test case
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% H formulation/s-formulation using g_top
i = any(g_top.faces.neighbors==0, 2);  % find all outer faces
I = i(g_top.cells.faces(:,1));         % vector of all faces of all cells, true if outer
j = false(4,1);                     % mask, cells can at most have 6 faces,
%j(1:2)=true;  %   extract east, west, north, south
j(bc_side)=true;%   extract east
J = j(g_top.cells.faces(:,2));         % vector of faces per cell, true if E,W,N,S
bc_faces = g_top.cells.faces(I & J, 1);
 %depth * density * gravity   
bc_h = addBC([], bc_faces, 'pressure', g_top.faces.z(bc_faces)*1000*norm(gravity));
bc_h = rmfield(bc_h,'sat');
bc_h.h=zeros(size(bc_h.face));  

bc_s = addBC([], bc_faces, 'pressure', g_top.faces.z(bc_faces)*1000*norm(gravity), 'sat', [0 1]);  

%% Coupled formulation
i = any(g.faces.neighbors==0, 2);  % find all outer faces
I = i(g.cells.faces(:,1));         % vector of all faces of all cells, true if outer
j = false(6,1);                     % mask, cells can at most have 6 faces,
%j(1:2)=true;  %   extract east, west, north, south
j(bc_side)=true;%   extract east
J = j(g.cells.faces(:,2));         % vector of faces per cell, true if E,W,N,S
bc_faces = g.cells.faces(I & J, 1);
                              
bc_coupled = addBC([], bc_faces, 'pressure', g.faces.centroids(bc_faces,3)*1000*norm(gravity), 'sat', [0 1]);

%% 3D grid
i = any(g3D.faces.neighbors==0, 2);  % find all outer faces
I = i(g3D.cells.faces(:,1));         % vector of all faces of all cells, true if outer
j = false(6,1);                     % mask, cells can at most have 6 faces,
%j(1:2)=true;  %   extract east, west, north, south
j(bc_side)=true;%   extract east
J = j(g3D.cells.faces(:,2));         % vector of faces per cell, true if E,W,N,S
bc_faces = g3D.cells.faces(I & J, 1);
bc_3D = addBC([], bc_faces, 'pressure', (g3D.faces.centroids(bc_faces,3)+H/2)*1000*norm(gravity), 'sat', [0 1]);


%% Make wells
ind_well=double(g_top.cells.indexMap(sub2ind(g_top.cartDims,ind(1),max(ind(2),1))));
W_h = addWell([], g_top, rock2d, ind_well,'Type', 'rate','Val',rate, ...
               'Radius',0.1,'comp_i',[1 0]);
W_s = [];


W_s = addWell(W_s, g_top, rock2d, ind_well,'Type', 'rate','Val',rate, ...
             'Radius',0.1,'comp_i',[1 0]);
          
          
W_coupled = verticalWell([], g, rockCoupled, ind(1), ind(2), k_ix, 'Type', 'rate','Val',rate, ...
             'Radius',0.1,'comp_i',[1 0]);          
          
W_3D = verticalWell([], g3D, rock, ind(1), ind(2), k_ix,'Type', 'rate','Val',rate,...
               'Radius',0.1,'comp_i',[1 0]);


src_s = [];
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


%% Define solvers to solve the system with different discretisations
%

% 1) original version in h formulation
%%{
n=1;
problem{n}.g = g_top;
problem{n}.fluid      = initVEFluidHForm(g_top, 'mu', mu, 'rho', rho,'sr',sr,'sw',sw);
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
problem{n}.g = g_top;
problem{n}.fluid = initSimpleVEFluidSForm('mu' , mu , 'rho', rho, ...
                           'height'  , g_top.cells.H,...
                           'sr', [sr, sw]);                        
problem{n}.sol = initResSolVE(g_top, 0);
S = computeMimeticIP(g_top, rock2d); %T=computeTrans(g_top,rock2d);
problem{n}.S=S;
problem{n}.bc=bc_s;
problem{n}.src=src_s;

problem{n}.psolver =@(sol,fluid,W,bc, src) ...
solveIncompFlow(sol, g_top, S, fluid,'wells',W,'bc',bc, 'src', src);
%     solveIncompFlowVE_s(sol, g_top, S, fluid,'wells',W,'bc',bc, 'src', src);
   
problem{n}.W= W_s; %W_3D;
problem{n}.tsolver =@(sol,fluid,dt,W,bc, src)...
         implicitTransport(sol, g_top, dt, rock2d, fluid,'vert_avrg',true, ...
                             'wells',W,'bc',bc, 'src', src);
% implicitTransportVE(sol, g_top, dt, rock2d, fluid,'vert_avrg',true, ...
%                              'wells',W,'bc',bc, 'src', src);
problem{n}.compute_h=true;
problem{n}.compute_sat=false;
problem{n}.col='gs';


n=3;
%  3) TPFA and coupled 3D grid and not integration in vertical
%     direction


%problem{n} = problem{n-1};
problem{n}.g = g;
problem{n}.sol = initResSolVE(g, 0);
problem{n}.fluid = initCoupledVEFluid('mu' , mu , 'rho', rho, ...
                           'height'  , g.cells.H,...
                           'sr', [sr, sw], 'region3D', region3D, 'n', [1, 1], 'g', g);  
T = computeTrans(g, rockCoupled);   
problem{n}.bc=bc_coupled;
problem{n}.W=W_coupled;
problem{n}.psolver =@(sol,fluid,W,bc, src)...
                     incompTPFAVE_coupled(sol, g, T, fluid,'wells', W, ...
                                    'bc',bc, 'src', src, 'region3D', region3D);
problem{n}.tsolver =@(sol,fluid,dt,W,bc, src)...
implicitTransport(sol, g, dt, rockCoupled, fluid,'vert_avrg',true, ...
                            'wells',W,'bc',bc, 'src', src, 'region3D', region3D, 'verbose', true);
problem{n}.compute_h=true;
problem{n}.compute_sat=false;
problem{n}.col='b*';
problem{n}.src=src_s;



n=4;
%  4) TPFA with 3D grid 
problem{n}.g = g3D;
problem{n}.fluid = initCoreyFluid('mu' , mu , 'rho', rho, ...
                                 'sr', [sr, sw], 'kwm', [1 1], 'n', [1 1]);  
                              %not entirely correct
problem{n}.fluid.sat2height = @(state)[ sat2height(state.s, g_top, rock)];                        
                            
problem{n}.sol        = initResSol(g3D, 0);
T = computeTrans(g3D, rock);    
problem{n}.psolver =@(sol,fluid,W,bc, src)...
                     incompTPFA(sol, g3D, T, fluid,'wells', W, ...
                                    'bc',bc, 'src', src);
problem{n}.tsolver =@(sol,fluid,dt,W,bc, src)...
implicitTransport(sol, g3D, dt, rock, fluid, ...
                            'wells',W,'bc',bc, 'src', src, 'verbose', true);
                             
problem{n}.bc=bc_3D;
problem{n}.src=src_s;

problem{n}.W=W_3D;
problem{n}.compute_h=true;
problem{n}.compute_sat=false;
problem{n}.col='c--';                                


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% initialize all problems with init solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for kk=1:numel(problem)
   problem{kk}.sol.s=zeros(problem{kk}.g.cells.num,1); %s_init;
   problem{kk}.sol.h=zeros(g_top.cells.num,1);
   problem{kk}.sol.s_max=problem{kk}.sol.s; %s_init;
   problem{kk}.sol.max_h=problem{kk}.sol.h; %h_init;
   problem{kk}.p_time=0;
   problem{kk}.t_time=0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run transport-simulation:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3);
t = 0;
while t < total_time
   figure(3),clf
   hold on 
   for kk=1:numel(problem)
     disp(kk) 
        
            free_volume= sum(problem{kk}.sol.h.*rock2d.poro.*g_top.cells.volumes);
      %end
      disp(['Total volume free ', num2str(free_volume/1e6)]);                    
    
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
      problem{kk}.sol = problem{kk}.psolver(problem{kk}.sol,problem{kk}.fluid,W,bc, src);
      tmp=toc(tmp);
      problem{kk}.p_time=problem{kk}.p_time+tmp;
      %% transport solve
      tmp=tic;
      nn=1;    
       for ii=1:nn         
         problem{kk}.sol = problem{kk}.tsolver(problem{kk}.sol,problem{kk}.fluid,dt/nn,W,bc, src);      
       end      
       
      tmp=toc(tmp);
      problem{kk}.t_time=problem{kk}.t_time+tmp;
      
      if(problem{kk}.compute_h)
         %[h,h_max]=problem{kk}.fluid.sat2height(problem{kk}.sol);
          [h]=problem{kk}.fluid.sat2height(problem{kk}.sol);
          if kk == 3            
             h = sat2height(problem{kk}.sol.s, g_top_coupled, rock);
             problem{kk}.sol.h = h;
             s = height2Sat(problem{kk}.sol, g_top, problem{kk}.fluid);
             problem{kk}.sol.s3D = s; 
             
             problem{kk}.sol.s3D(g.cells.inx3D(region3D)) =  problem{kk}.sol.s(region3D); 
          end
          
         problem{kk}.sol.h = h;
         problem{kk}.sol.max_h = h; %_max;
      end
      if(problem{kk}.compute_sat)         
         problem{kk}.sol.s =  (problem{kk}.sol.h*(1-sw)+(problem{kk}.sol.max_h-problem{kk}.sol.h)*sr)./g_top.cells.H;
         problem{kk}.sol.s_max =  problem{kk}.sol.max_h*(1-sw)./g_top.cells.H;
      end      
         
     
         subplot(2, 2, 1)
         hold on 
         plot(g_top.cells.centroids(:,1), problem{kk}.sol.h,problem{kk}.col)
         if(kk==1);axh=axis();else axis(axh);end
         ylim([0 H])
         title('Height')
         
         
         if kk == 2
            subplot(2, 2, 2)
            s = height2Sat(problem{kk}.sol, g_top, problem{kk}.fluid);
            plotCellData(g3D, s)
            view([0 0])
            caxis([0 1]); axis tight
            title(['years', num2str(t/year)])
            xlabel('Mimetic sform');
            if(g3D.cartDims(2)>1), view(3), end
         end
         
         if kk == 3            
        
         subplot(2, 2, 3)
         plotCellData(g3D, problem{kk}.sol.s3D)
         view([0 0])
         caxis([0 1]); axis tight
         xlabel('TPFA coupled 3D');
          if(g3D.cartDims(2)>1), view(3), end
         elseif kk == 4
           
         subplot(2, 2, 4)
         plotCellData(problem{kk}.g, problem{kk}.sol.s)
         view([0 0])
          caxis([0 1]), axis tight
          xlabel('TPFA full 3D');
          
         if(g3D.cartDims(2)>1), view(3), end
         
         end        
         
     
         
   end
   t=t+dt;
end

% s-formulation
figure
s = height2Sat(problem{2}.sol, g_top, problem{2}.fluid);
plotCellData(g3D, s)
view([0 0])
caxis([0 1]); axis tight off;
% title(['years', num2str(t/year)])
if(g3D.cartDims(2)>1), view(3), end

if printFigs
   set(gcf, 'render', 'painters')
   shading flat;
   print -depsc2 fig-form1D.eps
end

% coupled formulation
figure
% plotCellData(problem{kk}.g, problem{kk}.sol.s)
plotCellData(g3D, problem{3}.sol.s3D)
view([0 0])
caxis([0 1]); axis tight off;
if(g3D.cartDims(2)>1), view(3), end

if printFigs
   set(gcf, 'render', 'painters')
   shading flat;
   print -depsc2 fig-coup1D.eps
end


% 3D
figure
plotCellData(problem{4}.g, problem{kk}.sol.s)
view([0 0])
caxis([0 1]), axis tight off;
if(g3D.cartDims(2)>1), view(3), end

if printFigs
   set(gcf, 'render', 'painters')
   shading flat;
   print -depsc2 fig-org1D.eps
end
