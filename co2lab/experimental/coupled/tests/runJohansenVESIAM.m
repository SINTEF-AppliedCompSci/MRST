% Run coupled model 3D/VE for injection period, then run VE for post injection, 
% runs VE with s-formulation
% 

clear, clc; %close all;
disp('================================================================');
disp('   Vertical averaging applied to the Johansen formation');
disp('================================================================');
disp('');
try
   disp(' -> Reading JohansenVEgridSIAM.mat');
   load JohansenVEgridSIAM;
catch me
   disp(' -> Reading failed, constructing grid models');
   makeJohansenVEgridSIAM;
end



%% Set time and fluid parameters
coupled3D = true;
T          = 40*year();
stopInject = 40*year();
dT         = year()/5;
dTplot     = 1*dT;

% Make fluid structures
% at p = 300 bar
muw = 0.30860;  rhow = 975.86; %sw    = 0.1; 
muc = 0.056641; rhoc = 686.54; %srco2 = 0.2;
%kwm = [0.2142 0.85];
kwm = [1 1];
sw = 0;
srco2 = 0;


 fluidVE = initSimpleVEFluid('mu' , [muc muw] .* centi*poise, ...
                             'rho', [rhoc rhow] .* kilogram/meter^3, ...
                             'sr', [srco2, sw], 'height', G_top.cells.H);                         
                          
% NB: endre?                          
% fluid = initCoreyFluid('mu' , [muc muw] .* centi*poise, ...
%                        'rho', [rhoc rhow] .* kilogram/meter^3, ...
%                         'sr', [srco2 sw], 'kwm', kwm, 'n', [1 1]);                          
fluid = initCoupledVEFluid('mu' , [muc muw] .* centi*poise, ...
                       'rho', [rhoc rhow] .* kilogram/meter^3, ...
                        'sr', [srco2 sw], 'region3D', region3D, 'n', [1, 1], 'g', G_c);    
                     
fluid_eq = initCoupledVEFluid('mu' , [muw muw] .* centi*poise, ...
                       'rho', [rhow rhow] .* kilogram/meter^3, ...
                        'sr', [srco2 sw], 'region3D', region3D, 'n', [1, 1], 'g', G_c);                        
                     
                     
gravity on
                       
%% Set well and boundary conditions
% We use one well placed in the center of the model, perforated in layer 6.
% Injection rate is 1.4e4 m^3/day of supercritical CO2. Hydrostatic
% boundary conditions are specified on all outer boundaries that are not in
% contact with the shales; the latter are assumed to be no-flow boundaries.

% Set well in 3D model
%wellIx = [51, 51, 6, 6];

wellIx = [48, 51, 6, 6];


rate = 1.4e4*meter^3/day;
W = verticalWell([], G, rock, wellIx(1), wellIx(2), wellIx(3):wellIx(4),...
   'Type', 'rate', 'Val', rate, 'Radius', 0.1, 'comp_i', [1,0], 'name', 'I');


bc_c = addBC([], bcIx_c, 'pressure', G_c.faces.centroids(bcIx_c, 3)*rhow*norm(gravity), ...
            'sat', [0, 1]);
bc_c_nograv = addBC([], bcIx_c, 'pressure', G_c.faces.centroids(bcIx_c, 3)*rhow*norm(gravity)*0, ...
            'sat', [0, 1]);



% Well and BC in 2D model 
wellIxVE = find(G_top.columns.cells == W(1).cells(1));
wellIxVE = find(wellIxVE - G_top.cells.columnPos >= 0, 1, 'last' );
WVE      = addWell([], G_top, rock2D, wellIxVE, ...
                  'Type', 'rate','Val',rate,'Radius',0.1, 'comp_i',[1 0]);

                             
bcVE = addBC([], bcIxVE, 'pressure', G_top.faces.z(bcIxVE)*rhow*norm(gravity), 'sat', [0 1]);                  
                  

%% Prepare simulations
% Compute inner products and instantiate solution structure
%SVE = computeMimeticIPVE(G_top, rock2D, 'Innerproduct','ip_simple');
SVE = computeMimeticIP(G_top, rock2D, 'Innerproduct','ip_simple');
preComp = initTransportVE(G_top, rock2D); 
sol = initResSol(G_top, 0);
sol.wellSol = initWellSol(W, 300*barsa());
sol.h = zeros(G_top.cells.num, 1);
sol.max_h = sol.h; 
sol.s_max = sol.s;
%sol.s = height2Sat(sol, G_top, fluidVE);
%sol.s_max = sol.s; %height2Sat(sol, G_top, fluidVE);
%% Prepare plotting
% We will make a composite plot that consists of several parts:
%
% * a 3D plot of the plume
% * a pie chart of trapped versus free volume
% * a plane view of the plume from above
% * two cross-sections in the x/y directions through the well

% Extract grid to represent the 2D cross-sections
clear ijk
[ijk{1:3}] = ind2sub(G.cartDims, G.cells.indexMap);
ijk = [ijk{:}];

cells_sub_x = find(ijk(:,2) == wellIx(2));
G_x = extractSubgrid(G, cells_sub_x);
G_x = computeGeometry(G_x); 
G_x.cartDims = G.cartDims;

cells_sub_y = find(ijk(:,1) == wellIx(1));
G_y = extractSubgrid(G, cells_sub_y);
G_y = computeGeometry(G_y); 
G_y.cartDims = G.cartDims;
if ~coupled3D

% Figure for VE plot
scrsz  = get(0,'ScreenSize');
%figVE = figure('Position',[scrsz(3)-1024 scrsz(4)-700 1024 700]);
%set(0,'CurrentFigure',figVE);

% Saturation plotted in a 3D grid with the x- and y-slices outlined in cyan
% and magenta colours, respectively.
subplot(2,3,1:2); 
  plotGrid(G, 'FaceColor', 'none', 'EdgeAlpha', 0.1);
  plotGrid(G_x, 'FaceColor', 'none', 'EdgeColor', 'c');
  plotGrid(G_y, 'FaceColor', 'none', 'EdgeColor', 'm');
  plotWell(G, W, 'height', 500, 'color', 'r');
  title(['2D: CO2 saturation at ', ...
     num2str(convertTo(0,year)), ' years']);
  axis tight off, view([-85, 70]);
  caxis([0 1]);

% Pie chart of trapped and free CO2
subplot(2,3,3);
  pie([eps 1],{'trapped','free'});
  title('Total volume: 0');
  
% Height of the CO2 plume plotted on the 2D top-surface grid.
subplot(2,3,4);
  plotGrid(G_top, 'FaceColor', 'none', 'EdgeAlpha', 0.1);
  axis tight off, set(gca,'zdir','reverse'); view([-90, 90])
  title('Height of CO2-column');
  
% Saturation along the x-slice
subplot(4,3,8:9);
  plotGrid(G_x,'FaceColor','none','EdgeAlpha',0.1);
  axis tight off, view([0 0])
  
% Saturation along the y-slice
subplot(4,3,11:12);
  plotGrid(G_y,'FaceColor','none','EdgeAlpha',0.1);
  axis tight off, view([90 0])
end
%% Main loop
% Run the simulation using a sequential splitting with pressure and
% transport computed in separate steps. The transport solver is formulated
% with the height of the CO2 plume as the primary unknown and the relative
% height (or saturation) must therefore be reconstructed. 
[hsVE1,hsVE2] = deal([]);
t = 0;
fprintf(1,'\nSimulating %d years on injection',convertTo(stopInject,year));
fprintf(1,' and %d years of migration\n', convertTo(T-stopInject,year));
fprintf(1,'Time: %4d years', convertTo(t,year));

 cellNo  = rldecode(1:G_c.cells.num, diff(G_c.cells.facePos), 2) .';
plot_var  = @(x) plotCellData(G_c, x);
plot_div = @(x) plot_var(accumarray(cellNo, ...
   faceFlux2cellFlux(G_c, x.flux)));
   
if coupled3D
tic   

g_comptop = topSurfaceGrid(G_c);
T = computeTrans(G_c, rock_c);   

   sol3D = initResSol(G_c, 0);
   sol3D.wellSol = initWellSol(W, 300*barsa());
   sol3D.s_max = sol3D.s; 
   S = computeMimeticIP(G_c, rock_c, 'Innerproduct','ip_simple');
gravity on

soltmp = sol3D;
while t<stopInject
   % Advance solution: compute pressure and then transport
   
   gravity on
   sol3D = incompTPFAVE_coupled(sol3D, G_c, T, fluid,'wells', W, ...
       'bc',bc_c, 'region3D', region3D, 'verbose', true);
    % figure; plotFaces(G_c, G_c.facesBnd.index, sol3D.flux( G_c.facesBnd.index))
   %{
    gravity on
     soltmp = incompTPFAVE_coupled(soltmp, G_c, T, fluid_eq,'wells', W, ...
       'bc',bc_c, 'region3D', region3D, 'verbose', true);
     figure; plotFaces(G_c, G_c.facesBnd.index, soltmp.flux( G_c.facesBnd.index))
    
    %}
    %figure(5)
    %clf
    %plot_div(sol3D)
   
    gravity on
 sol3D = implicitTransportVE_coupled(sol3D, G_c, dT, rock_c, fluid,'vert_avrg',true, ...
          'wells',W,'bc',bc_c, 'region3D', region3D, 'verbose', true, 'nltol', 1e-5);
%   sol3D = implicitTransportVE_coupled(sol3D, G_c, dT, rock_c, fluid,'vert_avrg',true, ...
%           'wells',W,'bc',bc_c, 'region3D', region3D, 'verbose', true, 'nltol', 1e-6);
%     sol3D = implicitTransport(sol3D, G_c, dT, rock_c, fluid, ...
%       'bc', bc_c, 'wells', W);
%    
   % Reconstruct 'saturation' on real 3D grid    
  sol3D.h = sat2height(sol3D.s, g_comptop, rock_c);
  sol3D.s3D = height2Sat(sol3D, g_comptop, fluid);
  % sol.s(G_c.cells.inx3D(G_c.cells.region3D)) = sol3D.s(G_c.cells.region3D);
   % Check if we are to stop injecting
   if t>= stopInject
      WVE  = []; bcVE = []; 
   end
   fprintf(1,'Time: %4d years', convertTo(t,year));
   assert( max(sol.s(:,1))<1+eps && min(sol.s(:,1))>-eps );
   
   figure(9)
   clf
  subplot(2, 1, 1)
     plotGrid(G_c, 'FaceColor', 'none', 'EdgeAlpha', 0.1);
     plotCellData(G_c, sol3D.s, find(sol3D.s>0.0001), 'EdgeAlpha', 0.1);
     plotGrid(G_c, G_c.facesBnd.cells2D);
   
     title(['2D: CO2 saturation at ', ...
        num2str(convertTo(t, year)), ' years']);
  subplot(2, 1, 2)    
     plotGrid(g_comptop, 'FaceColor', 'none', 'EdgeAlpha', 0.1);
     %plotGrid(G_c, G_c.facesBnd.cells2D);
     plotCellData(g_comptop, sol3D.h,  find(sol3D.h(:,1)>5e-4), 'EdgeAlpha', 0.1);
          
   drawnow

   
   
   t = t + dT;
   
   % Plotting
   
end

%if mod(t,dTplot)~= 0, continue; end
   

%assert(all(sol3D.s(~G_c.cells.region3D)<1e-5))

% initialize VE simulation
%sol.h = sat2height(sol.s, G_top, rock);





%% Make figures
figure;  

plotGrid(g_comptop, 'FaceColor', 'none', 'EdgeAlpha', 0.1);
plotGrid(g_comptop, G_c.cells.mapTopSurface(G_c.facesBnd.cells2D), 'faceColor', 'none', 'edgeColor', 'r')
     %plotGrid(G_c, G_c.facesBnd.cells2D);
     plotCellData(g_comptop, sol3D.h,  find(sol3D.h(:,1)>5e-4), 'EdgeAlpha', 0.1);
     axis tight off
          view([-93 90])
          
          
          
figure;  

plotGrid(G, 'FaceColor', 'none', 'EdgeAlpha', 0.1);
plotFaces(G_c, G_c.facesBnd.index, 'faceColor', 'none', 'edgeColor', 'r', 'LineWidth',1.5)
%plotGrid(G, G_c.cells.inx3D(G_c.facesBnd.cells2D), 'faceColor', 'none', 'edgeColor', 'r', 'LineWidth',1.5)

cells3D = G_c.cells.inx3D(g_comptop.columns.cells(g_comptop.cells.columnPos(1:end-1)));

sol.h = zeros(G.cells.num,1);
sol.h(cells3D) = sol3D.h;

colorbar

     %plotGrid(G_c, G_c.facesBnd.cells2D);
     plotCellData(G, sol.h, find(sol.h(:,1)>5e-4), 'EdgeAlpha', 0.1);
     axis tight off
          view([-94 68])         
  
%if printFigs
%set(gcf, 'render', 'painters')
% shading flat; 
%  print -dpng fig-johansen40years_3D_2.png
%end


copySol = sol.h;

toc
end
freeVol = 0;
trappedVol = 0;
totVol = 0;


while t<T
   % Advance solution: compute pressure and then transport
   sol = solveIncompFlowVE_s(sol, G_top, SVE, fluidVE, ...
      'bc', bcVE, 'wells', WVE);
%    sol = explicitTransportVE(sol, G_top, dT, rock, fluidVE, ...
%       'bc', bcVE, 'wells', WVE, 'preComp', preComp,'dt_dynamic',true);
    sol = implicitTransportVE(sol, G_top, dT, rock2D, fluidVE, ...
      'bc', bcVE, 'wells', WVE);
   
   % Reconstruct 'saturation' defined as s=h/H, where h is the height of
   % the CO2 plume and H is the total height of the formation
   sol.h = fluid.sat2height(sol);
   sol.s3D = height2Sat(sol, G_top, fluidVE);
   assert( max(sol.s(:,1))<1+eps && min(sol.s(:,1))>-eps );
   t = t + dT;
   
   % Check if we are to stop injecting
   if t>= stopInject
      WVE  = []; bcVE = []; 
   end
   
   % Compute trapped and free volumes of CO2
   fprintf(1,'\b\b\b\b\b\b\b\b\b\b%4d years', convertTo(t,year));
%    freeVol = ...
%       sum(sol.h.*rock2D.poro.*G_top.cells.volumes)*(1-fluidVE.sw);
%    trappedVol = ...
%       sum((sol.max_h-sol.h).*rock2D.poro.*G_top.cells.volumes)*fluidVE.sr;
%    totVol = trappedVol + freeVol;

   % Plotting
   if mod(t,dTplot)~= 0, continue; end
   set(0,'CurrentFigure',figVE);
   subplot(2,3,1:2);
     delete(hsVE1);
     hsVE1 = plotCellData(G, sol.s3D, find(sol.s3D>0.01), 'EdgeAlpha', 0.1);
     title(['2D: CO2 saturation at ', ...
        num2str(convertTo(t, year)), ' years']);
     caxis([0 1]);
%    subplot(2,3,3);
%      str1 = ['trapped ' num2str(round(trappedVol/totVol*100)) ' %'];
%      str2 = ['free ' num2str(round(freeVol/totVol*100)) ' %'];
%      pie([max(trappedVol,eps) freeVol],{str1, str2}); 
%      title(['Total volume: ', ...
%         num2str(round(convertTo(totVol,mega*meter^3))),' M m^3']);
   subplot(2,3,4);
     delete(hsVE2);
     hsVE2 = plotCellData(G_top, sol.h, ...
        find(sol.h(:,1)>5e-3),'EdgeAlpha', 0.1);
   subplot(2,3,5); cla
     plotCellData(G_x, sol.s3D(cells_sub_x),'EdgeColor','k','EdgeAlpha',0.1);
     axis tight off, view([0 0])
     title('CO2-saturation, x-slice');
   subplot(2,3,6); cla
     plotCellData(G_y, sol.s3D(cells_sub_y),'EdgeColor','k','EdgeAlpha',0.1);
     axis tight off, view([90 0])
     title('CO2-saturation, y-slice');
   drawnow
end
fprintf(1,'\n\n');

displayEndOfDemoMessage(mfilename)

%save VEsim.mat sol3D copySol sol