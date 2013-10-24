clear;  close all


load JohansenVEgrid;
g = G; g_top = G_top; bcIx_v = bcIxVE;
clear G G_top bcIxVE;

 
gravity on
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Params to change:
% options for vertical average transport solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

int_vert= true; 
int_vert_poro= false; %true; %false;
semi_implicit=false; %true; %true;

% time:                          
numYears = 510;                 
T      = numYears*year();
dT     = T/(numYears);
dTplot = 2*dT; %100*dT; %10*year();
stopInject = 110*year;
run3D  = false; %true;
plotting = true;
only2D = false;
dire = 'x';

% well
wellIx = [51, 51, 6, 6];
rate = 1.4e4*meter^3/day;

%% %%%%%%% SUBGRID FOR SIM/PLOT %%%%%%%%%%%%%%%%%%%%%%
clear ijk
[ijk{1:3}] = ind2sub(g.cartDims, g.cells.indexMap);
ijk = [ijk{:}];

cells_sub_x = find(ijk(:,2) == wellIx(2));
g_sub_x = extractSubgrid(g, cells_sub_x);
g_sub_x = computeGeometry(g_sub_x); g_sub_x.cartDims = g.cartDims;
cells_sub_y = find(ijk(:,1) == wellIx(1));
g_sub_y = extractSubgrid(g, cells_sub_y);
g_sub_y = computeGeometry(g_sub_y); g_sub_y.cartDims = g.cartDims;
 
if only2D  % Simulate in 2D
  nx = g.cartDims(1); ny=g.cartDims(2); nz=g.cartDims(3); 
  nx_v = g_top.cartDims(1); ny_v=g_top.cartDims(2); nz_v=1;
  if strcmp(dire, 'y')
     g = g_sub_y; g_top = topSurfaceGrid(g);
     rock.perm = rock.perm(cells_sub_y,:);
     rock.poro = rock.poro(cells_sub_y);
   %  wellIx = [1 wellIx(2) wellIx(3) wellIx(4)];    
     ix1 = boundaryFaceIndices(g, 'BACK', 1:nx, 1:nz, 1:ny);
     %ix2 = boundaryFaceIndices(g, 'FRONT', 1:nx, 1:nz, 1:ny);
     ix1_v = boundaryFaceIndices(g_top, 'BACK', 1:nx, 1, []);
     %ix2_v = boundaryFaceIndices(g_top, 'FRONT', 1:nx, 1, []);   
  else
     g = g_sub_x; g_top = topSurfaceGrid(g);
     rock.perm = rock.perm(cells_sub_x,:);
     rock.poro = rock.poro(cells_sub_x);
    % wellIx = [wellIx(1) 1 wellIx(3) wellIx(4) ];    
     %ix1 = boundaryFaceIndices(g, 'LEFT', 1:ny, 1:nz, 1:nx);     
     ix1 = boundaryFaceIndices(g, 'RIGHT', 1:ny, 1:nz ,  1:nx);
     %ix1_v = boundaryFaceIndices(g_top, 'LEFT', 1:ny,  1, []);
     ix1_v = boundaryFaceIndices(g_top, 'RIGHT', 1:ny,  1, []);     
  end
  bcIx = [ix1]; % ix2];  
  bcIx_v = [ix1_v]; % ix2_v];  
  
  rate = rate/10;
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make fluid structures

% at p = 300 bar
muw = 0.30860; muc = 0.056641; rhow = 975.86; rhoc = 686.54;
kwm = [0.2142 0.85];
srco2 = 0.2; sw = 0.1;

% 3D:
fluid = initCoreyFluid('mu' , [muc muw] .* centi*poise, ...
                        'rho', [rhoc rhow] .* kilogram/meter^3, ...
                        'n'  , [2, 2],  'sr', [srco2 sw], 'kwm', kwm);
% 2D:
fluid_v = initVEFluid(g_top, 'mu' , [muc muw] .* centi*poise, ...
                           'rho', [rhoc rhow] .* kilogram/meter^3, ...
                          'sr', srco2, 'sw', sw, 'kwm', kwm);
%% %%%%%%%%%%% Make averaged version of rock %%%%%%%%
rock2D  = averageRock(rock, g_top);
%rock2D.perm = rock.perm(g_top.cells.map3D,1);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set well and boundary data - hydrostatic boundary conditions on all sides
% 3D 
W = verticalWell([], g, rock, wellIx(1), wellIx(2), wellIx(3):wellIx(4),  'Type', 'rate', 'Val', rate,  ...
      'Radius', 0.1, 'comp_i', [1,0], 'name', 'I');
pressure = sum(g.faces.centroids(bcIx,:).*repmat(gravity(), length(bcIx),1) ,2)*rhow;                
bc = addBC([], bcIx, 'pressure', pressure, 'sat', [0 1]);

% 2D 
ind_well1 = find(g_top.columns.cells == W(1).cells(1));
ind_well1 = find(ind_well1-g_top.cells.columnPos>=0, 1, 'last' );
W_v = addWell([], g_top, rock2D, ind_well1, 'Type', 'rate','Val',rate,'Radius',0.1);

bc_v = addBC([], bcIx_v, 'pressure', g_top.faces.z(bcIx_v)*rhow*norm(gravity));
bc_v = rmfield(bc_v,'sat');
bc_v.h=zeros(size(bc_v.face));

% for 2D/vertical average, we need to change defintion of the wells
for i=1:numel(W_v)
   W_v(i).compi=nan;
   W_v(i).h = g_top.cells.H(W_v(i).cells);
end


%% Make inner products 
% 3D:
if run3D
   S  = computeMimeticIP(g, rock, 'Verbose', true);
end
% 2D/vertical average
S_2d = computeMimeticIPVE(g_top, rock2D,'Innerproduct','ip_simple');
% precompute transport quantities
preComp = initTransportVE(g_top, rock2D); 

%% Make initial solution structure
% 3D:
rSol         = initResSol (g, 350*barsa, 0.0); rSol.wellSol = initWellSol(W, 300*barsa());
% 2D
sol = initResSol(g_top, 0); sol.wellSol = initWellSol(W, 300*barsa());
sol.h = zeros(g_top.cells.num, 1); sol.max_h = sol.h; 

pv     = poreVolume(g,rock);
pv_3D(g_top.columns.cells)=rock.poro(g_top.columns.cells)...
         .*rldecode(g_top.cells.volumes,diff(g_top.cells.columnPos)); %*(1-fluid_v.sw);

      
      
if plotting % Prepare plotting of saturations   
   if~only2D
      if run3D
         figure(12)
         plotGrid(g, 'FaceColor', 'none', 'EdgeAlpha', 0.1);
         plotWell(g, W, 'height', 500, 'color', 'w');
         axis off, view(3) %,
         title(['3D: CO2 saturation at ', num2str(convertTo(0,year)), ' years']);  caxis([0 1])
      else
         figure(12)
         plotGrid(g_top, 'FaceColor', 'none', 'EdgeAlpha', 0.1);
         axis off, view(2)
      end
      figure(13)
      plotGrid(g, 'FaceColor', 'none', 'EdgeAlpha', 0.1);
      plotWell(g, W, 'height', 500, 'color', 'w');
      title(['2D: CO2 saturation at ', ...
         num2str(convertTo(0,year)), ' years']);
      axis off, view(3);
      caxis([0 1]);  
   else
      
      if strcmp(dire, 'y'), v = [90 13]; else v = [0 0];   end
      figure(15) %subplot(2, 1, 1)
      plotGrid(g, 'FaceColor', 'none', 'EdgeAlpha', 0.1); drawnow; view(v);axis tight;
      figure(16); % subplot(2, 1, 2)
      plotGrid(g, 'FaceColor', 'none', 'EdgeAlpha', 0.1); drawnow; view(v);axis tight;
      drawnow; view(v) ;
      title(['CO2 saturation at ', num2str(convertTo(0,year)), ' years']);
   end
end

hs2 = []; ha2=[];  hs = []; ha=[];
% Start the main loop
t  = 0;  plotNo = 1;

bhp2D = []; bhp3D = []; error = []; dt_3D = 0; dt_mrst = []; dt_va = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Main loop alternate between solving the transport and the flow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while t < T,
   if run3D
      % 3D: ------------------------------------------------------------------
      % Solve pressure equation
      rSol = solveIncompFlow(rSol, g, S, fluid, 'wells', W, 'bc', bc);
      % Solve transport equation
      [rSol,report] = explicitTransportOld(rSol, g, dT, rock, fluid, 'wells', W, 'bc', bc);
      % Check for inconsistent saturations
      assert(max(rSol.s(:,1)) < 1+eps && min(rSol.s(:,1)) > -eps);
      
      dt_mrst = [dt_mrst; dT/report.timesteps];
   end
   
   % 2D: ------------------------------------------------------------------
   % Solve pressure equation 
   
   sol = solveIncompFlowVE(sol, g_top, S_2d, rock, fluid_v, ...
      'bc', bc_v,'wells',W_v);
   
   % Solve transport equation
   if run3D, vararg = {'dt', dt_3D}; else vararg = {}; end
   [sol, dt_2D] = explicitTransportVE(sol, g_top, dT, rock, fluid_v,...
      'Verbose',false,'intVert',int_vert, 'intVert_poro', int_vert_poro, ...
      'semi_implicit', semi_implicit,'bc', bc_v,'wells',W_v,'preComp', preComp); %, vararg{:});
   
   dt_va = [dt_va; dt_2D];
   % convert from height to saturation
   sol.s = height2Sat(sol, g_top, fluid_v);
   %-----------------------------------------------------------------------
   bhp3D =   [bhp3D; rSol.wellSol.pressure];
   bhp2D = [bhp2D; sol.wellSol.pressure];
   assert(max(sol.s(:,1)) < 1+eps && min(sol.s(:,1)) > -eps);
   % Increase time:
   t = t + dT;
   
   error = [error abs(sum(sol.s.*g.cells.volumes.*rock.poro)-rate*t)/(rate*t)];
      %sum(rSol.s(:,1).*g.cells.volumes.*rock.poro)
   
     
   %%%%%%%%%%%%%%%%%%%%%% Fra testSleipnerPaulo %%%%%%%%%%%%%%%%%%%%%%%%%%%
   if  mod(t, dTplot) == 0           
      disp(['Time  ',num2str(t/year)])
      if(int_vert_poro)
         free_volume = sum(integrateVertically(pv_3D', sol.h, g_top))*(1-fluid_v.sw);       
      else
         free_volume= sum(sol.h.*rock2D.poro.*g_top.cells.volumes)*(1-fluid_v.sw);
      end
      disp(['Total volume free ', num2str(free_volume/1e6)]);
      if(isfield(sol,'max_h'))
         if(int_vert_poro)
            trapped_volume= sum(integrateVertically(pv_3D', sol.max_h, g_top)- ...
                                integrateVertically(pv_3D', sol.h, g_top))*fluid_v.sr;
         else
            trapped_volume=sum((sol.max_h-sol.h).*rock2D.poro.*g_top.cells.volumes)*fluid_v.sr;
         end
         disp(['Total volume trapped ', num2str(trapped_volume/1e6)]);
         disp(['Total volume ', num2str((trapped_volume+free_volume)/1e6)]);
         %disp(['Total volume 3D: 'num2str(rate*t/1e6)])
         totrate = min(rate*t, rate*stopInject);
         disp(['Mass consv h: ',num2str((trapped_volume+free_volume-totrate)/(totrate))])        
         disp(['Mass consv sat:', num2str((sum(sol.s.*pv)-totrate)/(totrate))])
        % disp(['Mass consv sat:', num2str((sum(sol.s.*pv))/1e6)])
         disp('******************')
      end     
   end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   if t >= stopInject && (~isempty(bc_v) ||~isempty(W_v)) %stop injection
      disp('Stop injection')
      W_v = [];   bc_v = []; W = [];  bc = []; %should perhaps keep bc?
   %  
   end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   if plotting
      if mod(t, dTplot) == 0           
         if ~only2D
         if run3D
            figure(12); delete(hs)
            %clf; hs= plotCellData(g, rSol.pressure); cx = caxis;
            hs= plotCellData(g, rSol.s(:,1),find(rSol.s(:,1)>0.01), 'EdgeAlpha', 0.1);
            view(3),  caxis([0 1]); drawnow; title(['3D: CO2 saturation at ', num2str(convertTo(t,year)), ' years']);
         else
            figure(12); delete(hs)
           hs = plotCellData(g_top, sol.h(:,1),find(sol.h(:,1)>0.01), 'EdgeAlpha', 0.1);
           caxis([0 max(g_top.cells.H)])
            colorbar
         end
         figure(13); delete(hs2)
         %clf; hs2 = plotCellData(g_top, sol.pressure); caxis(cx)
         hs2 = plotCellData(g, sol.s,find(sol.s>0.01), 'EdgeAlpha', 0.1);
         title(['2D: CO2 saturation at ', num2str(convertTo(t,year)), ' years']);
         caxis([0 1]); view(3); drawnow         
         figure(15) % 2D plots
         subplot(2, 2, 1)
         title(['3D: CO2 saturation at ', num2str(convertTo(t,year)), ' years']);
         plotCelldata(g_sub_x, rSol.s(cells_sub_x,1)); caxis([0 1]); drawnow; view([0 0]);axis tight;
         subplot(2, 2, 2)
         title(['2D: CO2 saturation at ', num2str(convertTo(t,year)), ' years']);
         plotCelldata(g_sub_x, sol.s(cells_sub_x)); caxis([0 1]); drawnow; view([0 0]);axis tight
          subplot(2, 2, 3)
         plotCelldata(g_sub_y, rSol.s(cells_sub_y,1)); caxis([0 1]); drawnow; view([90 0]);axis tight;
          subplot(2, 2, 4)
         plotCelldata(g_sub_y, sol.s(cells_sub_y)); caxis([0 1]); drawnow; view([90 0]);axis tight;
         else                
            if strcmp(dire, 'y'), v = [90 13]; else v = [0 0];   end
           if run3D, figure(15);  delete(hs)          
             hs =  plotCelldata(g, rSol.s, find(rSol.s(:,1)>0.01),'EdgeAlpha', 0.1 ); caxis([0 1]); drawnow; view(v);axis tight;end
           figure(16);   delete(hs2)
             hs2 =  plotCelldata(g,sol.s, find(sol.s>0.01),'EdgeAlpha', 0.1 ); caxis([0 1]); drawnow; view(v) ;   axis tight;    
                title(['CO2 saturation at ', num2str(convertTo(t,year)), ' years']);
         end    
      end
   end         
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


if ~plotting
   if mod(t, dTplot) == 0      
      if ~only2D
         if run3D
            figure(12); %delete(hs)
            plotGrid(g, 'FaceColor', 'none', 'EdgeAlpha', 0.1);
            %clf; hs= plotCellData(g, rSol.pressure); cx = caxis;
            hs= plotCellData(g, rSol.s(:,1),find(rSol.s(:,1)>0.01), 'EdgeAlpha', 0.1);
            view(3),  caxis([0 1]); drawnow; title(['3D: CO2 saturation at ', num2str(convertTo(t,year)), ' years']);
         else
            figure(12); %delete(hs)
            plotGrid(g, 'FaceColor', 'none', 'EdgeAlpha', 0.1);
            plotCellData(g_top, sol.h(:,1),find(sol.h(:,1)>0), 'EdgeAlpha', 0.1)
            caxis([0 max(g_top.cells.H)])
            colorbar
         end
         figure(13); %delete(hs2)
         plotGrid(g, 'FaceColor', 'none', 'EdgeAlpha', 0.1);
         %clf; hs2 = plotCellData(g_top, sol.pressure); caxis(cx)
         hs2 = plotCellData(g, sol.s,find(sol.s>0.01), 'EdgeAlpha', 0.1);
         title(['2D: CO2 saturation at ', num2str(convertTo(t,year)), ' years']);
         caxis([0 1]); view(3); drawnow
         
         figure(15) % 2D plots
         subplot(2, 2, 1)
         plotCelldata(g_sub_x, log10(rock.perm(cells_sub_x,1))); drawnow; view([0 0]);axis tight;
         subplot(2, 2, 3)
         plotCelldata(g_sub_y, log10(rock.perm(cells_sub_y)));  drawnow; view([90 0]);axis tight;
         subplot(2, 2, 2)
         plotCelldata(g_sub_x, sol.s(cells_sub_x)); caxis([0 1]); drawnow; view([0 0]);axis tight;
         subplot(2, 2, 4)
         plotCelldata(g_sub_y, sol.s(cells_sub_y)); caxis([0 1]); drawnow; view([90 0]);axis tight;
      else         
         if strcmp(dire, 'y'), v = [90 13]; else v = [0 0];   end
         figure(15);  delete(hs)
         plotGrid(g, 'FaceColor', 'none', 'EdgeAlpha', 0.1);
         hs =  plotCelldata(g, rSol.s, find(rSol5.s(:,1)>0.01),'EdgeAlpha', 0.1 ); caxis([0 1]); drawnow; view(v);axis tight;
         figure(16);   delete(hs2)
         plotGrid(g, 'FaceColor', 'none', 'EdgeAlpha', 0.1);
         hs2 =  plotCelldata(g,sol.s, find(sol.s>0.01),'EdgeAlpha', 0.1 ); caxis([0 1]); drawnow; view(v) ;   axis tight;
         title(['CO2 saturation at ', num2str(convertTo(t,year)), ' years']);
      end
   end
end


disp(['CFL 3D max:', num2str(max(dt_mrst/year))])
disp(['CFL 3D min:', num2str(min(dt_mrst/year))])
disp(['CFL VA max:', num2str(max(dt_va/year))])
disp(['CFL VA min:', num2str(min(dt_va/year))])


%disp(error')
% if plotting
%    figure;
%    plotCellData(g, sol.s,find(sol.s>0.01), 'EdgeAlpha', 0.1);
%    plotWell(g, W, 'height', 500, 'color', 'w');
% end
