clc; clear; close all;
%% Preliminaries
gravity reset
gravity on
T          = 500*year();
stopInject = 100*year();
dT         = .01*year();%.1*year();
dTplot     = 1*dT;
injectTop = false;

rate = 1.4e4*meter^3/day;

% Well radius
r = 0.1;
% Fluid data at p = 300 bar
muw = 0.30860;  rhow = 975.86; sw    = 0.1;
muc = 0.056641; rhoc = 686.54; srco2 = 0.2;
kwm = [0.2142 0.85];

v = [-30 20];
dims = [30, 1, 3];
% dims = [10, 10, 10];
mrstModule add co2lab-matlab/vertical-equil-lab
%% Create grids
% Create a logically Cartesian grid and add in geometry data
grdecl = simpleGrdecl(dims, 0, 'undisturbed', true, 'physDims', [1000 1000 10]);
G_full = processGRDECL(grdecl);
%%
clear ijk;
[ijk{1:3}] = ind2sub(G_full.cartDims, G_full.cells.indexMap(:));
ijk        = [ijk{:}];
%%
mi = max(ijk(:,1));
% % Define a 3D region in the middle of the domain
r3D = ismember(ijk(:,1), round(mi/3):round(2*mi/3) );
% Define a 2D region in the middle of the domain
% r3D = ~ismember(ijk(:,1), round(mi/3):round(2*mi/3) );
clf;
plotGrid(G_full, 'facea', 0)
plotGrid(G_full, r3D)

% G = cartGrid(dims, dims*100);

% Perturb the z coordinates based on the x coordinates
perturb = @(x) .025*(max(x) - min(x))*sin(x)  + .25*(x-min(x));


G_full.nodes.coords(:,3) = G_full.nodes.coords(:,3); + perturb(G_full.nodes.coords(:,1));

G_full = computeGeometry(G_full);
% Create top surface grid
[G_top, G_full] = topSurfaceGrid(G_full);



figure(11),clf;
subplot(3,1,1)
plotGrid(G_full);
subplot(3,1,2)
plotGrid(G_top);
%% Show grids
h = figure(1); clf;
% Show the full 3D grid
subplot(2,1,1)
plotGrid(G_full);
title('Full grid')
axis equal tight; view(v)
% set(gca, 'ZDir', 'normal')
% Show the top grid
subplot(2,1,2)
plotGrid(G_top, 'FaceColor', 'red')
title('Surface grid')
axis equal tight; view(v)
% set(gca, 'ZDir', 'normal')
%% Fluid and petrophysical properties
rock.perm = 25*milli*darcy*ones(G_full.cells.num,1);
% rock.perm(:,2) = rock.perm(:,2)./1000;
% rock.perm = (1 + 24*rand(G.cells.num,1)).*milli*darcy;
rock.poro = 0.3*ones(G_full.cells.num,1);

rock2D    =    averageRock(rock, G_top);





% S_2D       = computeMimeticIPVE(G_2D, rock2D, 'Innerproduct','ip_simple');
S_2D       = computeMimeticIP(G_top, rock2D, 'Innerproduct','ip_simple');
%S          = computeMimeticIP(G_full, rock, 'Innerproduct','ip_simple');
T_full          = computeTrans(G_full, rock);
T_2D = computeTrans(G_top, rock2D);
preComp = initTransportVE(G_top, rock2D);

%%
clear ijk
[ijk{1:3}] = ind2sub(G_full.cartDims, G_full.cells.indexMap(:));
y_middle = double(median(double(ijk{2})));
x_min    = double(min(ijk{1}));
x_max    = double(max(ijk{1}));
% Injector in the left side of the domain
wI = ijk{1} == x_min + 3 & ijk{2} == y_middle;
if injectTop
    wI = wI & ijk{3} == max(ijk{3});
end
% Producer in the right side of the domain
wP       = ijk{1} == x_max - 3 & ijk{2} == y_middle;
% Injector
% W = addWell([], G_full, rock, find(wI),...
%    'Type', 'rate', 'Val', rate, 'Radius', r, 'comp_i', [1,0], 'name', 'I');
% % Producer
% W = addWell(W, G_full, rock, find(wP),...
%    'Type', 'bhp', 'Val', 0, 'Radius', r, 'comp_i', [1,0], 'name', 'P');

W = verticalWell([], G_full, rock, x_min, y_middle, [],...
   'Type', 'rate', 'Val', rate, 'Radius', r, 'comp_i', [1,0], 'name', 'I');
% Producer
W = verticalWell(W, G_full, rock, x_max, y_middle, [],...
   'Type', 'bhp', 'Val', 0, 'Radius', r, 'comp_i', [1,0], 'name', 'P');



% Convert the wells to 2D equivialent
W_2D = convertwellsVE(W, G_full, G_top, rock2D);




clf; 
plotGrid(G_full, 'edgea', .1, 'facea', .1);
% plotWell(G, W, 'radius', 10);
plotGrid(G_full, vertcat(W.cells), 'facec', 'red');
axis equal tight; view(v)
%%
mu = [muc muw] .* centi*poise; rho = [rhoc rhow] .* kilogram/meter^3; 
fluid_H = initVEFluidHForm(G_top, 'mu' , mu, ...
                                 'rho', rho, ...
                                 'sr', srco2, 'sw', sw, 'kwm', kwm);

fluid_S = initSimpleVEFluidSForm('mu' , mu , 'rho', rho, ...
                                'height'  , G_top.cells.H,...
                                'sr', [srco2, sw]); 
                                              
fluid_full = initCoreyFluid('mu' , mu , 'rho', rho, 'sr', [srco2, sw], 'kwm', [1 1], 'n', [1 1]);  
%%
% H formulation
sol_H = initResSolVE(G_top, 0);
sol_H.wellSol = initWellSol(W, 0);
% S formulation
sol_S = initResSolVE(G_top, 0);
sol_S.wellSol = initWellSol(W, 0);

% Full grid
sol_full = initResSol(G_full, 0);
sol_full.wellSol = initWellSol(W, 0);


% % %% Prepare plot
% % 
% opts = {'slice', double([x_min + 3 y_middle]), 'Saxis', [0 1-fluid_H.sw], 'maxH', max(G.nodes.coords(:,3)), ...
%         'view', v, 'wireH', true, 'wireS', true};
% % % plotPanelVE(G, G_2D, W, sol_S, 0.0, [0 0 1], opts{:});
% plotPanelVE(G, G_2D, W, sol_H, 0.0, [0 0 1], opts{:});
%% Main loop
% Run the simulation using a sequential splitting with pressure and
% transport computed in separate steps. 
t = 0;
fprintf(1,'\nSimulating %d years on injection',convertTo(stopInject,year));
fprintf(1,' and %d years of migration\n', convertTo(T-stopInject,year));
fprintf(1,'Time: %4d years', convertTo(t,year));
tic;

while t<T
    
   % Advance solution:
   % First we compute pressures
   sol_H =    solveIncompFlowVE  (sol_H, G_top, S_2D, rock, fluid_H, 'wells', W_2D);
   sol_S =    solveIncompFlowVE_s(sol_S, G_top, S_2D,       fluid_S, 'wells', W_2D);
   sol_S.pressure = sol_S.pressure./G_top.cells.H;
   %sol_full = solveIncompFlow(sol_full , G_full, S, fluid_full, 'wells', W);
   sol_full = incompTPFA(sol_full , G_full, T_full, fluid_full, 'wells', W);
   
   
   % which is then used to advance saturations
   sol_H = explicitTransportVE(sol_H, G_top, dT, rock, fluid_H, ...
                              'wells', W_2D, ...
                              'preComp', preComp, 'Verbose', false, 'computeDt', true);
                          
   
   sol_S       = implicitTransport(sol_S, G_top, dT, rock2D, fluid_S, 'wells', W_2D);

   sol_full    = implicitTransport(sol_full, G_full, dT, rock, fluid_full, 'wells', W);
   
   
   % Reconstruct 'saturation' defined as s=h/H, where h is the height of
   % the CO2 plume and H is the total height of the formation
   sol_H.s = height2Sat(sol_H, G_top, fluid_H);
   % For the s formulation, construct height from the saturation values
   sol_S.h = fluid_S.sat2height(sol_S);
   %... monotone?
   sol_S.h_max = max(sol_S.h, sol_S.h_max);
   
   %
   sol_full.h = sat2height(sol_full.s, G_top, rock);
   
   
   assert( max(sol_H.s(:,1))<1+eps && min(sol_H.s(:,1))>-eps );
   
   
   
   
   
   t = t + dT;

   % Check if we are to stop injecting
   if t>= stopInject
      W_2D  = []; bcVE = []; dT = 5*year(); dTplot = dT;
   end

   % Compute trapped and free volumes of CO2
   fprintf(1,'\b\b\b\b\b\b\b\b\b\b%4d years', int32(convertTo(t,year)));
%    fprintf(1, '%d years\n', int32(convertTo(t, year)));
   vol = volumesVE(G_top, sol_H, rock2D, fluid_H);
   
%    figure(1); clf; hold on
%    x = G_top.cells.centroids(:,1);
%    plot(x, [sol_S.h sol_H.h, sol_full.h])
%    legend({'S formulation', 'H formulation', 'Full 3d'})
%    axis([0 max(G_full.nodes.coords(:,1)) 0 max(G_top.cells.H)])
%    drawnow
   figure(2);
   subplot(2,2,1)
   
   [s h] = normalizeValuesVE(G_top, sol_S, fluid_S);
   plotCellData(G_full, s)
%    plotCellData(G_full, height2Sat(sol_S, G_top, fluid_S))
   axis square tight; view([0 0])
   title('S formulation')  
   
   subplot(2,2,2)
   [s h] = normalizeValuesVE(G_top, sol_H, fluid_H);
   plotCellData(G_full, s)
   plotCellData(G_full, height2Sat(sol_H, G_top, fluid_H))
   axis square tight; view([0 0])
   title('H formulation')
   
   subplot(2,2,3)
   plotCellData(G_full, sol_full.s)
   title('Full 3D')
   axis square tight; view([0 0])
   
   drawnow
   %%
   a=figure(22),clf,hold on
   legend={};
   plot(G_top.cells.centroids(:,1),G_top.cells.z,'r','LineWidth',2)
   plot(G_top.cells.centroids(:,1),G_top.cells.z+G_top.cells.H,'r','LineWidth',2)
   [s h] = normalizeValuesVE(G_top, sol_S, fluid_S)
   plot(G_top.cells.centroids(:,1),G_top.cells.z+h)
   %plot(G_top.cells.centroids(:,1),G_top.cells.z+sol_S.s.*G_top.cells.H)
   legend{end+1}={'S formulationVE'}
   [s h] = normalizeValuesVE(G_top, sol_H, fluid_H);
   plot(G_top.cells.centroids(:,1),G_top.cells.z+h)
   legend{end+1}={'H formulationVE'}
   set(gca,'YDir','reverse')
  %%
   
   figure(3);
   subplot(2,2,1)
   plotCellData(G_top, sol_S.pressure)
   axis square tight; view(0,90)
   title('S formulation')
   
   subplot(2,2,2)
   plotCellData(G_top, sol_H.pressure)
   axis square tight; view(0,90)
   title('H formulation')
   
   subplot(2,2,3)
   plotCellData(G_full, sol_full.pressure)
   title('Full 3D')
   axis square tight; view(0,90)

   
   
   drawnow;
   
   

end
fprintf(1,'\n\n');
%% Check fluxes etc
bndf = G_coupled.cells.faces(G_coupled.facesBnd.cellFace3D);
% sol_coupled.flux(bndf)
sum(sol_coupled.flux(bndf).*G_coupled.faces.areas(bndf))
mean(sol_coupled.flux(sol_coupled.flux>0))

% cellNo = rldecode(1:G_coupled.cells.num, double(diff(G_coupled.cells.facePos)), 2).';

figure(5); clf;
% f = boundaryFaces(G);
% f = find(G_coupled.faces.tag > 0);
f = setdiff(find(G_coupled.faces.normals(:,3) == 0), boundaryFaces(G_coupled) );
flux = sol_coupled.flux./G_coupled.faces.areas;
% flux = T_coupled;
% flux = flux(G_coupled.cells.faces(:,1));
plotFaces(G_coupled, f, flux(f), 'facea', .2);
colorbar;
std(flux(f))

figure(6); clf;
% f = boundaryFaces(G);
% f = find(G_coupled.faces.tag > 0);
flux = sol_full.flux;
f = setdiff(find(G_full.faces.normals(:,3) == 0), boundaryFaces(G_full) );
% flux = computeTrans(G_full, rock);
% flux = flux(G_full.cells.faces(:,1));
plotFaces(G_full, f, flux(f), 'facea', .2);
colorbar;
std(flux(f))
%%
tarea = T_coupled./G_coupled.faces.areas(G_coupled.cells.faces(:,1));
tt = unique(tarea);
figure(7); clf;

plotFaces(G_coupled, G_coupled.cells.faces((tarea == tt(1)), 1), 'FaceC', 'r')
% plotFaces(G_coupled, G_coupled.cells.faces((tarea == tt(2)), 1), 'FaceC', 'g')
% plotFaces(G_coupled, G_coupled.cells.faces((tarea == tt(3)), 1), 'FaceC', 'b')
