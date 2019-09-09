%% run field model
clear;close all;clc
mrstModule add incomp mimetic  ad-blackoil ad-core glpk ntpfa_glpk ...
    ad-props mpfa eni mrst-gui nfvm ntpfa blackoil-sequential ...
    wellpaths
plotconc = false;

% myfile = 'grid_data/tet2.dat';
% G = importGrid(myfile,'tetBench');
% G = twister(G, 0.1);
% G = computeGeometry(G);

% nx = 7;
% ny = 5;
% nz = 3;
% G = cartGrid([nx, ny, nz]);
% G = twister(G, 0.11);
% G = computeGeometry(G);

nx = 7;
ny = 5;
nz = 3;
G = cartGrid([nx, ny, nz], [1 1 1]);
%G = tetrahedralGrid(G.nodes.coords);
%G = twister(G, 0.1);
G = computeGeometry(G);

figure
plotGrid(G,'facealpha', 0.5)

% Extract faces in x direction
tol = 1e-5;
pts = [tol, 0.5+tol, 0.5+tol;
       1-tol, 0.5+tol, 0.5+tol];
wellpath = makeSingleWellpath(pts);
%plotWellPath(wellpath);
cells_xdir = findWellPathCells(G, wellpath);
[Gsub,~,subgf] = extractSubgrid(G, cells_xdir);

rock = makeRock(G, 100*milli*darcy, 1);
pv = sum(poreVolume(G,rock));

mu_value = 1;
rho_value = 1;
fluid = initSingleFluid('mu',mu_value,'rho',rho_value);
fluidad = initSimpleADIFluid('phases', 'W', 'mu' , mu_value, 'rho', ...
                             rho_value);
fluid_glpk = initSimpleADIFluid();
p0 = 15e6; % from FlowNTPFA call below
s0 = 1.0;

% Wells
Wtp = [];
cellsWell1 = well_cells(G, '1');
cellsWell2 = well_cells(G, '2');

radius = 0.05/10;
% Wtp = addWell(Wtp, G, rock, cellsWell1, 'Type', 'rate', ...
%             'InnerProduct', 'ip_tpf', ...
%             'comp_i', [1], ...
%             'Val', 1.0/day(), 'Radius', radius, 'Name', 'I');
% Wtp = addWell(Wtp, G, rock, cellsWell2, 'Type', 'bhp', ...
%             'InnerProduct', 'ip_tpf', ...
%             'comp_i', [1], ...
%             'Val', 1.0e5, 'Radius', radius, 'Dir', 'y', 'name', 'P');

% Wtp = addWell(Wtp, G, rock, cellsWell1, 'Type', 'bhp', ...
%             'InnerProduct', 'ip_tpf', ...
%             'comp_i', [1], ...
%             'Val', -1.0e5, 'Radius', radius, 'Dir', 'y', 'name', 'P');

% bc_nfvm.face=boundaryFaces(G);
% bc_nfvm.type=repmat({'flux'},[numel(bc_nfvm.face),1]);
% bc_nfvm.value=repmat({@(x)0},[numel(bc_nfvm.face),1]);
% bc_nfvm.face=boundaryFaces(G);
% bc_nfvm.type=repmat({'pressure'},[numel(bc_nfvm.face),1]);
% bc_nfvm.value=repmat({@(x)0},[numel(bc_nfvm.face),1]);

% bc_std = [];
% bc_std = pside(bc_std, G, 'xmin', 1);
% bc_std = pside(bc_std, G, 'xmax', 1);
% bc_std = pside(bc_std, G, 'ymin', 1);
% bc_std = pside(bc_std, G, 'ymax', 1);
% bc_std = pside(bc_std, G, 'zmin', 1);
% bc_std = pside(bc_std, G, 'zmax', 1);

% % Must hard code these?
neumann_val = 1; % Hard coded in convertBC
dirichlet_val = 0; % Hard coded in convertBC

% Don't use fluxside / pside but set manually to get correct
% scaling. 
bf = boundaryFaces(G);
bc_std = [];
isbdry = (1:G.faces.num)';
tol = 1e-5;
xmin = min(G.nodes.coords(:,1));
xmax = max(G.nodes.coords(:,1));
imin = G.faces.centroids(:,1) < xmin+tol;
imax = G.faces.centroids(:,1) > xmax-tol;
neumann_ix = isbdry(imin);
area = G.faces.areas(neumann_ix);
bc_std = addBC(bc_std, neumann_ix, 'flux', neumann_val.*area, 'sat', []);
dirichlet_ix = isbdry(imax);
bc_std = addBC(bc_std, dirichlet_ix, 'pressure', dirichlet_val, ...
               'sat', []);
% Set the rest as hom Neumann
bc_std = addBC(bc_std, setdiff(bf, union(neumann_ix, dirichlet_ix)), ...
               'flux', 0.0, 'sat', []);


% Convert the traditional BC struct to useful format for ntpfa codes
bc_flow_ntpfa = convertBC2FlowNTPFA(G, bc_std);
bc_nfvm_new = bc_std;
bc_glpk = bc_std;
bc_glpk.sat = [1, 0];

%state0 = initState(G, Wtp, p0);
%state0 = initResSol(G, p0, 1);
% mfd:
state0 = initState(G, Wtp, p0, s0);
state0_glpk = initResSol(G, p0, [0 s0]);

% Schedule with dummy dt
dt = 1;
schedule = simpleSchedule(dt, 'W', Wtp, 'bc', bc_std);

% Set saturation manually
%bc_std.sat = s0*ones(numel(bc_std.face), 1);
%bc_glpk.sat = [ones(numel(bc_glpk.face), 1), ...
%               zeros(numel(bc_glpk.face), 1)];

solvers = {};
results = {};



% %% MPFA Xavier style
% solvers{end+1} = 'mpfa';
% disp(solvers{end});
% %mpfastruct = computeNeumannMultiPointTrans(G, rock);
% %state = incompMPFA3(G, mpfastruct, Wtp, 'outputFlux', true);
% mpfastruct = computeMultiPointTrans2(G, rock);
% state = incompMPFAbc(G, mpfastruct, bc_std, 'outputFlux', true);
% figure,plotCellData(G,state.pressure/1e6);plotWell(G,Wtp);view(3)
% title(['Pressure ', solvers{end}]);h=colorbar;h.Label.String='Pressure [Pa]';clear h;
% if plotconc
%     Qinj=sum(state.wellSol(1).flux);tt=pv/Qinj;nstep=100;dt=tt/nstep;
%     state = tracerTransport_implicit(G,rock,Wtp,state,dt,nstep);
%     figure, plotCellData(G,state.cres(:,end));plotWell(G,Wtp);
%     view(3);colorbar;title(['Concentration ', solvers{end}]);
%     drawnow
% end
% results{end+1} = state;


% %% TPFA from NFVM
% solvers{end+1} = 'tpfa';
% disp(solvers{end});
% state=FlowTPFA(G,TransTPFA(G,rock,bc_flow_ntpfa),fluid,Wtp);
% figure,plotCellData(G,state.pressure/1e6);plotWell(G,Wtp);view(3)
% title(['Pressure ', solvers{end}]);h=colorbar;h.Label.String='Pressure [Pa]';clear h;
% Qinj=sum(state.wellSol(1).flux);tt=pv/Qinj;nstep=100;dt=tt/nstep;
% state=tracerTransport_implicit(G,rock,Wtp,state,dt,nstep);
% figure, plotCellData(G,state.cres(:,end));plotWell(G,Wtp);
% view(3);colorbar;title(['Concentration ', solvers{end}]);
% drawnow
% results{end+1} = state;

%% TPFA mrst version
solvers{end+1} = 'mrst tpfa';
disp(solvers{end});
T = computeTrans(G,rock);
state = incompTPFA(state0, G, T, fluid, 'wells', Wtp, 'bc', bc_std);
figure,plotToolbar(G,state);plotWell(G,Wtp);view(3)
title(['Pressure ', solvers{end}]);h=colorbar;h.Label.String='Pressure [Pa]';clear h;
figure,plotFaceData(G,cells_xdir, state.flux);colorbar,view(3)
title(['flux ', solvers{end}])
if plotconc
    Qinj=sum(state.wellSol(1).flux);tt=pv/Qinj;nstep=100;dt=tt/nstep;
    state = tracerTransport_implicit(G,rock,Wtp,state,dt,nstep);
    figure, plotCellData(G,state.cres(:,end));plotWell(G,Wtp);
    view(3);colorbar;title(['Concentration ', solvers{end}]);
drawnow
end
results{end+1} = state;


%% nonlinear TPFA
solvers{end+1} = 'NFVM Kobaisi';
disp(solvers{end});
interpFace=findHAP(G,rock,bc_flow_ntpfa);
disp(['fraction of faces with centroids outside convex hull ', num2str(interpFace.percentage)]);
interpFace=correctHAP(G,interpFace);
OSflux=findOSflux(G,rock,bc_flow_ntpfa,interpFace);
state=FlowNTPFA(G,bc_flow_ntpfa,fluid,Wtp,OSflux,p0*ones(G.cells.num,1),1e-7,1000);
%figure,plotCellData(G,state.pressure/1e6);view(3)
figure,plotToolbar(G,state);plotWell(G,Wtp);view(3)
title(['Pressure ', solvers{end}]);h=colorbar;h.Label.String='Pressure [Pa]';clear h;
figure,plotFaceData(G,cells_xdir, state.flux);colorbar,view(3)
title(['flux ', solvers{end}])
if plotconc
    Qinj=sum(state.wellSol(1).flux);tt=pv/Qinj;nstep=100;dt=tt/nstep;
    state=tracerTransport_implicit(G,rock,Wtp,state,dt,nstep);
    figure, plotCellData(G,state.cres(:,end));plotWell(G,Wtp);
    view(3);colorbar;title(['Concentration ', solvers{end}]);
    drawnow
end
results{end+1} = state;



% %% nonlinear NFVM in new framework
% solvers{end+1} = 'NFVM new';
% disp(solvers{end});
% model = GenericBlackOilModel(G, rock, fluidad,'water', true, 'oil', false, 'gas', false);
% %model = WaterModel(G, rock, fluidad);
% model = model.validateModel();
% model.FluxDiscretization.PermeabilityPotentialGradient.PermeabilityGradientDiscretization = NFVM(model, bc_nfvm_new);
% model.FacilityModel = model.FacilityModel.setupWells(Wtp);
% [wellSols, states] = simulateScheduleAD(state0, model, schedule);
% %figure,plotCellData(G,states{1}.pressure/1e6);plotWell(G,Wtp);view(3)
% figure,plotToolbar(G,states);plotWell(G,Wtp);view(3)
% title(['Pressure ', solvers{end}]);h=colorbar;h.Label.String='Pressure [Pa]';clear h;
% if plotconc
%     Qinj = sum(states{1}.wellSol(1).flux);tt=pv/Qinj;nstep=100;dt=tt/nstep;
%     state = tracerTransport_implicit(G,rock,Wtp,states{1},dt,nstep);
%     figure,plotCellData(G,state.cres(:,end));plotWell(G,Wtp);
%     view(3);colorbar;title(['Concentration ', solvers{end}]);
% end
% results{end+1} = states{1};


% % NTPFA glpk
% solvers{end+1} = 'NTPFA-glpk';
% disp(solvers{end});
% mrstVerbose on
% mrstDebug on
% model = PressureOilWaterModelNTPFAopt(G,rock,fluid_glpk);
% state = incompSinglePhaseNTPFA(model,state0_glpk,'bc', bc_glpk,'src', []);
% figure,plotToolbar(G,state);plotWell(G,Wtp);view(3)
% title(['Pressure ', solvers{end}]);h=colorbar;h.Label.String='Pressure [Pa]';clear h;
% if plotconc
%     Qinj=sum(state.wellSol(1).flux);tt=pv/Qinj;nstep=100;dt=tt/nstep;
%     state = tracerTransport_implicit(G,rock,Wtp,state,dt,nstep);
%     figure, plotCellData(G,state.cres(:,end));plotWell(G,Wtp);
%     view(3);colorbar;title(['Concentration ', solvers{end}]);
%     drawnow
% end
% results{end+1} = state;




% %% TPFA new framework
% solvers{end+1} = 'TPFA new';
% disp(solvers{end});
% model = GenericBlackOilModel(G, rock, fluidad,'water', true, 'oil', false, 'gas', false);
% model = model.validateModel();
% model.FluxDiscretization.PermeabilityPotentialGradient.PermeabilityGradientDiscretization ...
%     = TwoPointFluxApproximation(model);
% %model.FacilityModel = model.FacilityModel.setupWells(Wtp);
% [wellSols, states] = simulateScheduleAD(state0, model, schedule);
% %figure,plotCellData(G,states{1}.pressure/1e6);plotWell(G,Wtp);view(3)
% figure,plotToolbar(G,states);plotWell(G,Wtp);view(3)
% title(['Pressure ', solvers{end}]);h=colorbar;h.Label.String='Pressure [Pa]';clear h;
% Qinj=sum(states{1}.wellSol(1).flux);tt=pv/Qinj;nstep=100;dt=tt/nstep;
% state = tracerTransport_implicit(G,rock,Wtp,states{1},dt,nstep);
% figure, plotCellData(G,state.cres(:,end));plotWell(G,Wtp);
% view(3);colorbar;title(['Concentration ', solvers{end}]);
% results{end+1} = state;



% %% nonlinear MPFA (slow and not accurate?)
% solvers{end+1} = 'nonlinear MPFA';
% disp(solvers{end});
% %state=NewtonMPFA(G,fluid,Wtp,OSflux,1e-9,100);
% state=FlowNMPFA(G,bc_flow_ntpfa,fluid,Wtp,OSflux,1e-9,10000);
% figure,plotCellData(G,state.pressure/1e6);plotWell(G,Wtp);axis image;view(3)
% title('Pressure-NMPFA');h=colorbar;h.Label.String='Pressure [Pa]';clear h;
% Qinj=sum(state.wellSol(1).flux);tt=1.5*pv/Qinj;nstep=100;dt=tt/nstep;
% state=tracerTransport_implicit(G,rock,Wtp,state,dt,nstep);
% figure, plotCellData(G,state.cres(:,end));plotWell(G,Wtp);
% view(3);colorbar;title(['Concentration ', solvers{end}]);
% results{end+1} = state;


%% MFD
solvers{end+1} = 'mfd';
disp(solvers{end});
state=incompMimetic(state0,G,computeMimeticIP(G,rock),fluid,'wells',Wtp,'bc',bc_std);
%figure, plotCellData(G,state.pressure/1e6);view(3);
figure,plotToolbar(G,state);plotWell(G,Wtp);view(3)
title(['Pressure ', solvers{end}]);h=colorbar;h.Label.String='Pressure [Pa]';clear h;
figure,plotFaceData(G,cells_xdir, state.flux);colorbar,view(3)
title(['flux ', solvers{end}])
if plotconc
    Qinj=sum(state.wellSol(1).flux);tt=pv/Qinj;nstep=100;dt=tt/nstep;
    state=tracerTransport_implicit(G,rock,Wtp,state,dt,nstep);
    figure, plotCellData(G,state.cres(:,end));plotWell(G,Wtp);
    view(3);colorbar;title(['Concentration ', solvers{end}]);
    drawnow
end
results{end+1} = state;


% % primitive comparison of mfd and ntpfa_kobaisi
% figure
% plotCellData(G,abs(results{1}.pressure - results{2}.pressure)./results{2}.pressure); 
% colorbar 
% view(3)
% title('relative error pressure')
% disp(['max pressures ', num2str(max(results{1}.pressure),'%1.2e'), ' ', ...
%       num2str(max(results{2}.pressure), '%1.2e')]) 

% if plotconc
%     figure
%     plotFaceData(G,abs(results{1}.flux - results{2}.flux)./results{2}.flux);
%     colorbar 
%     view(3)
%     title('relative error flux')
%     disp(['max flux ', num2str(max(results{1}.flux),'%1.2e'), ' ', ...
%           num2str(max(results{2}.flux), '%1.2e')]) 
% end





%%
% data1=abs(state.pressure-state.pressure)./1e6;
% data2=abs(state.pressure-state.pressure)./1e6;
% xmin=min([data1;data2]);xmax=max([data1;data2]);
% figure,subplot(1,2,1);
% plotCellData(G,data1,'edgealpha',0.2);view(-80,50);axis tight off;
% h=colorbar('horiz');h.Label.String='Pressure [Pa]';h.Label.FontSize=15;
% title('\mid\itp_{\rmtpfa}-\itp_{\rmmfd}\rm\mid','fontsize',20);caxis([xmin xmax]);
% subplot(1,2,2),plotCellData(G,data2,'edgealpha',0.2);view(-80,50);axis tight off;
% h=colorbar('horiz');h.Label.String='Pressure [Pa]';h.Label.FontSize=15;
% title('\mid\it{p}_{\rmntpfa}-\it{p}_{\rmmfd}\rm\mid','fontsize',20);caxis([xmin xmax]);
% clear data1 data2 xmin xmax h
% 
% data=abs(state.cres(:,end)-state.cres(:,end));
% figure,subplot(1,2,1),
% plotCellData(G,data,'edgealpha',0.2);view(-80,50);axis tight off
% h=colorbar('horiz');h.Label.String='Concentration';h.Label.FontSize=15;caxis([0 1]);
% title('\mid\itC_{\rmtpfa}-\itC_{\rmmfd}\rm\mid','fontsize',20);
% data=abs(state.cres(:,end)-state.cres(:,end));
% subplot(1,2,2),plotCellData(G,data,'edgealpha',0.2);view(-80,50);axis tight off
% h=colorbar('horiz');h.Label.String='Concentration';h.Label.FontSize=15;caxis([0 1]);
% title('\mid\itC_{\rmntpfa}-\itC_{\rmmfd}\rm\mid','fontsize',20);clear data h

if plotconc
    figure, hold on
    for k = 1:numel(results)
        plot((1:nstep)'/nstep, ...
             results{k}.cwell, 'linewidth', 2);
    end
    legend(solvers)
    xlabel('pore volume injected');
    ylabel('concentration');
end
