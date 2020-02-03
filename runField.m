%% run field model
clear all;close all;clc
mrstModule add incomp mimetic  ad-blackoil ad-core glpk ntpfa_glpk ...
    ad-props mpfa project-eni mrst-gui ntpfa nfvm blackoil-sequential ...
    wellpaths 

plotconc = false;
plotjac = true;
plotflux = false;

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
refine = 1;
G = cartGrid([nx, ny, nz]*refine, [10 20 30]);
G = tetrahedralGrid(G.nodes.coords);
%G = twister(G, 0.1);
G = computeGeometry(G);

figure
plotGrid(G,'facealpha', 0.5)

% Extract faces in x direction for plotting fluxes
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
Wtp = addWell(Wtp, G, rock, cellsWell1, 'Type', 'rate', ...
            'InnerProduct', 'ip_tpf', ...
            'comp_i', [1], ...
            'Val', 1.0/day(), 'Radius', radius, 'Name', 'I');
Wtp = addWell(Wtp, G, rock, cellsWell2, 'Type', 'bhp', ...
            'InnerProduct', 'ip_tpf', ...
            'comp_i', [1], ...
            'Val', 1.0e5, 'Radius', radius, 'Dir', 'y', 'name', 'P');

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


% Case with Neumann on xmin and Dirichlet on xmax. Hom Neumann
% otherwise. 
%bc_std = bc_neumann_and_dirichlet(G);

% Case with only hom Neumann
bc_std = bc_hom_neumann(G);

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

% %% MPFA Xavier 
% solvers{end+1} = 'mpfa';
% disp(solvers{end});
% %mpfastruct = computeNeumannMultiPointTrans(G, rock);
% %state = incompMPFA3(G, mpfastruct, Wtp, 'outputFlux', true);
% mpfastruct = computeMultiPointTrans2(G, rock);
% state = incompMPFAbc(G, mpfastruct, bc_std, 'outputFlux', true);
% figure,plotCellData(G,state.pressure/1e6);plotWell(G,Wtp);view(3)
% title(['Pressure ', solvers{end}]);h=colorbar;h.Label.String='Pressure [Pa]';clear h;
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
state = incompTPFA(state0, G, T, fluid, 'wells', Wtp, 'bc', bc_std, ...
                   'matrixOutput', true);
figure,plotToolbar(G,state);plotWell(G,Wtp);view(3)
title(['Pressure ', solvers{end}]);h=colorbar;h.Label.String='Pressure [Pa]';clear h;
results{end+1} = state;


%% nonlinear TPFA
solvers{end+1} = 'NFVM Kobaisi';
mrstverbosestatus=mrstVerbose;
mrstVerbose on
disp(solvers{end});
interpFace=findHAP(G,rock,bc_flow_ntpfa);
disp(['fraction of faces with centroids outside convex hull ', num2str(interpFace.fraction)]);
interpFace=correctHAP(G,interpFace);
OSflux=findOSflux(G,rock,bc_flow_ntpfa,interpFace);
state=FlowNTPFA(G,bc_flow_ntpfa,fluid,Wtp,OSflux,p0*ones(G.cells.num,1),1e-7,1000,'matrixOutput', true);
%figure,plotCellData(G,state.pressure/1e6);view(3)
figure,plotToolbar(G,state);plotWell(G,Wtp);view(3)
title(['Pressure ', solvers{end}]);h=colorbar;h.Label.String='Pressure [Pa]';clear h;
results{end+1} = state;
mrstVerbose(mrstverbosestatus);


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


%% MFD (Must be at the end)
solvers{end+1} = 'mfd';
disp(solvers{end});
state=incompMimetic(state0,G,computeMimeticIP(G,rock),fluid,'wells',Wtp,'bc',bc_std,'matrixOutput', true);
%figure, plotCellData(G,state.pressure/1e6);view(3);
figure,plotToolbar(G,state);plotWell(G,Wtp);view(3)
title(['Pressure ', solvers{end}]);h=colorbar;h.Label.String='Pressure [Pa]';clear h;
figure,plotFaceData(G,cells_xdir, state.flux);colorbar,view(3)
title(['flux ', solvers{end}])
results{end+1} = state;


%% primitive comparisons with MFD (must be last)
L2err = @(uh,u,vol) num2str(sqrt(sum((uh-u).^2.*vol))); % ./ sqrt(sum(u.^2.*vol)));
for i = 1:numel(results)-1 % last is mfd
    disp(' ')
    disp(solvers{i})
    figure
    diff_pressure = abs(results{i}.pressure - results{end}.pressure)./abs(results{end}.pressure);
    plotCellData(G, diff_pressure); colorbar; view(3)
    titlename = ['relative error pressure ', solvers{i}, ' vs ', solvers{end}];
    title(titlename)
    disp(titlename)
    disp(['diff p max mean L2: ', ...
          num2str(max(diff_pressure)), ' ', ...
          num2str(mean(diff_pressure)), ' ', ...
          L2err(results{i}.pressure, results{end}.pressure, G.cells.volumes)])

    figure
    diff_flux = abs(results{i}.flux - results{end}.flux); %./abs(results{2}.flux);
    plotFaceData(G, diff_flux);colorbar; view(3)
    titlename = ['absolute error flux ', solvers{i}, ' vs ', solvers{end}];
    title(titlename)
    disp(titlename)
    disp(['diff flux max mean L2: ', ...
          num2str(max(diff_flux)), ' ', ...
          num2str(mean(diff_flux)), ' ', ...
          L2err(results{i}.flux, results{end}.flux, G.faces.areas)])
end



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

if plotflux
    for k = 1:numel(results)
        figure
        plotFaceData(G,cells_xdir, results{k}.flux);
        colorbar,view(3)
        title(['flux ', solvers{k}])
    end
end


if plotjac
    for k = 1:numel(results)
        if contains(solvers{k}, 'kobaisi', 'IgnoreCase', true)
            A = results{k}.jac;
            titlesuffix = ' Jacobian';
        else
            A = results{k}.A;
            titlesuffix = '';
        end

        figure
        spy(A)
        title([solvers{k}, titlesuffix])
    end
end


if plotconc
    for k = 1:numel(results)
        state = results{k};
        Qinj=sum(state.wellSol(1).flux);tt=pv/Qinj;nstep=100;dt=tt/nstep;
        state=tracerTransport_implicit(G,rock,Wtp,state,dt,nstep);
        figure, plotCellData(G,state.cres(:,end));plotWell(G,Wtp);
        view(3);colorbar;title(['Concentration ', solvers{k}]);
        drawnow
    end
    
    figure, hold on
    for k = 1:numel(results)
        plot((1:nstep)'/nstep, ...
             results{k}.cwell, 'linewidth', 2);
    end
    legend(solvers)
    xlabel('pore volume injected');
    ylabel('concentration');
end
