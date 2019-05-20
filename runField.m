% run field model
clear;close all;clc
mrstModule add incomp mimetic  ad-blackoil ad-core glpk ntpfa_glpk ...
    ad-props mpfa eni

% $$$ %% Grid and rock
% $$$ %grdecl=readGRDECL('NORNE.grdecl');
% $$$ grdecl = readGRDECL(fullfile(getDatasetPath('norne'), 'NORNE.GRDECL'));
% $$$ usys=getUnitSystem('METRIC');
% $$$ grdecl=convertInputUnits(grdecl,usys);
% $$$ G=processGRDECL(grdecl,'Tolerance',0.05);G=G(1);G=computeGeometry(G);
% $$$ rock=grdecl2Rock(grdecl,G.cells.indexMap);
% $$$ is_pos=rock.perm(:,3)>0;rock.perm(~is_pos,3)=1e-6*min(rock.perm(is_pos,3));
% $$$ is_pos=rock.perm(:,2)>0;rock.perm(~is_pos,2)=1e-6*min(rock.perm(is_pos,2));
% $$$ is_pos=rock.perm(:,1)>0;rock.perm(~is_pos,1)=1e-6*min(rock.perm(is_pos,1));
% $$$ figure, plotGrid(G);view(3);clear grdecl usys GG is_pos
% $$$ % ----------------------------------------------------------------------
% $$$ x=G.cells.centroids(:,1);y=G.cells.centroids(:,2);
% $$$ ind=x>4.59e5|y>7.324e6;[~,~,K]=ind2sub(G.cartDims,G.cells.indexMap);
% $$$ ind=ind|K<4;cells=find(ind);G=removeCells(G,cells);clear cells K
% $$$ rock.perm=rock.perm(~ind,:);rock.poro=rock.poro(~ind);rock.ntg=rock.ntg(~ind);
% $$$ % -------------------------------------------------------------------------------
% $$$ figure,subplot(1,2,1);
% $$$ plotCellData(G,log10(rock.perm(:,1)),'edgecolor','none');view(3);
% $$$ c=colorbar('southoutside');c.Label.String='log_{10}(\itk_{\rmH})';
% $$$ subplot(1,2,2);
% $$$ plotCellData(G,log10(rock.perm(:,3)),'edgecolor','none');view(3)
% $$$ c=colorbar('southoutside');c.Label.String='log_{10}(\itk_{\rmV})';clear c;
% $$$ % ----------------------------------------------------------------
% $$$ x=4.562e5;y=7.3205e6;
% $$$ d=bsxfun(@minus,G.cells.centroids(:,1:2),[x y]);
% $$$ d=sqrt(dot(d,d,2));[~,inj]=min(d);
% $$$ x=4.583e5;y=7.3234e6;
% $$$ d=bsxfun(@minus,G.cells.centroids(:,1:2),[x y]);
% $$$ d=sqrt(dot(d,d,2));[~,pro]=min(d);
% $$$ nz=G.cartDims(3);
% $$$ [I,J,~]=ind2sub(G.cartDims,G.cells.indexMap(inj));
% $$$ Wtp=verticalWell([],G,rock,I,J,1:nz,'type','bhp','val',20e6,'name','I','Comp_i', 1,'refDepth',min(G.cells.centroids(:,end)));
% $$$ Wm=verticalWell([],G,rock,I,J,1:nz,'type','bhp','val',20e6,'name','I','Comp_i', 1,'innerproduct','ip_quasirt','refDepth',min(G.cells.centroids(:,end)));
% $$$ [I,J,~]=ind2sub(G.cartDims,G.cells.indexMap(pro));
% $$$ Wtp=verticalWell(Wtp,G,rock,I,J,1:nz,'type','bhp','val',10e6,'name','P','Comp_i', 1,'refDepth',min(G.cells.centroids(:,end)));
% $$$ Wm=verticalWell(Wm,G,rock,I,J,1:nz,'type','bhp','val',10e6,'name','P','Comp_i', 1,'innerproduct','ip_quasirt','refDepth',min(G.cells.centroids(:,end)));
% $$$ clear x y d ind inj pro I J nz



% overwrite
nx = 10;
ny = 10;
nz = 5;
G = cartGrid([nx, ny, nz]);
G = computeGeometry(G);



    parentdir = 'data/two_faults_small/HYBRID_GRID/MATFILE';
    workdatadir = fullfile(parentdir, 'working_data');
    filename = fullfile(parentdir, 'HG_MRSTGrid.mat');
    load(filename);
    G = G_2;
    volume_hmin = min(G.cells.volumes.^(1/3));
    area_hmin = min(sqrt(G.faces.areas));
    hmin = min(volume_hmin, area_hmin);
    G = removePinch(G, max(eps, 1e-6*hmin));
    G = computeGeometry(G);

% $$$     % G = createHexTetGrid(nx, ny, nz, [nx, 0, 0]);
% $$$     G = tetrahedralGrid(G.nodes.coords);
% $$$     G = computeGeometry(G);

rock = makeRock(G, 100*milli*darcy, 1);


val1 = 20e6;
val2 = 10e6;
radius = 0.01;
cellsWell1 = well_cells(G, '1');
Wtp = addWell([], G, rock, cellsWell1, 'Type', 'bhp', ...
              'InnerProduct', 'ip_quasirt', ...
              'comp_i', [1], ...
              'Val', val1, 'Radius', radius, 'name', 'I');
cellsWell2 = well_cells(G, '2');
Wtp = addWell(Wtp, G, rock, cellsWell2, 'Type', 'bhp', ...
            'InnerProduct', 'ip_quasirt', ...
            'comp_i', [1], ...
            'Val', val2, 'Radius', radius, 'Dir', 'y', 'name', 'P');



pv=sum(poreVolume(G,rock));
bc.face=boundaryFaces(G);
bc.type=repmat({'flux'},[numel(bc.face),1]);
bc.value=repmat({@(x)0},[numel(bc.face),1]);

% load fieldModel;
figure, subplot(2,2,1),plotGrid(G);view(-5,55);plotWell(G,Wtp,'height',100);axis tight off
subplot(2,2,3),plotCellData(G,log10(convertTo(rock.perm(:,1),milli*darcy)),'edgecolor','none');
view(-5,55);h=colorbar;h.Label.String='log_{10}(\itk_{\rmH})';axis tight off;set(h,'fontsize',12)
subplot(2,2,4),plotCellData(G,log10(convertTo(rock.perm(:,end),milli*darcy)),'edgecolor','none');
view(-5,55);h=colorbar;h.Label.String='log_{10}(\itk_{\rmV})';axis tight off;set(h,'fontsize',12)
subplot(2,2,2),plotCellData(G,rock.poro,'edgecolor','none');
view(-5,55);h=colorbar;h.Label.String='\phi';axis tight off;set(h,'fontsize',12)
drawnow
mu_value = 1;
rho_value = 1;
fluid=initSingleFluid('mu',mu_value,'rho',rho_value);
%% linear TPFA
disp('tpfa')
stp=FlowTPFA(G,TransTPFA(G,rock,bc),fluid,Wtp);
figure,plotCellData(G,stp.pressure/1e6);plotWell(G,Wtp);view(3)
title('Pressure-TPFA');h=colorbar;h.Label.String='Pressure[MPa]';clear h;
Qinj=sum(stp.wellSol(1).flux);tt=pv/Qinj;nstep=100;dt=tt/nstep;
stp=tracerTransport_implicit(G,rock,Wtp,stp,dt,nstep);
figure, plotCellData(G,stp.cres(:,end));plotWell(G,Wtp);
view(3);colorbar;title('Concentration-TPFA');
drawnow

% mrst version
p0 = 15e6; % from FlowNTPFA call below
T = computeTrans(G,rock);
state0 = initState(G, Wtp, p0);
stpmrst = incompTPFA(state0, G, T, fluid, 'wells', Wtp);
stpmrst=tracerTransport_implicit(G,rock,Wtp,stpmrst,dt,nstep);

%% MFD
disp('mfd')
sm=incompMimetic(initState(G,Wtp,0),G,computeMimeticIP(G,rock),fluid,'wells',Wtp);
figure, plotCellData(G,sm.pressure/1e6);plotWell(G,Wtp);view(3);
title('Pressure-MFD');h=colorbar;h.Label.String='Pressure[MPa]';clear h;
Qinj=sum(sm.wellSol(1).flux);tt=pv/Qinj;nstep=100;dt=tt/nstep;
sm=tracerTransport_implicit(G,rock,Wtp,sm,dt,nstep);
figure, plotCellData(G,sm.cres(:,end));plotWell(G,Wtp);
view(3);colorbar;title('Concentration-MFD');
drawnow

%%
disp('nonlinear TPFA')
interpFace=findHAP(G,rock,bc);
disp(['faces with centroids outside convex hull ', num2str(interpFace.percentage)]);
interpFace=correctHAP(G,interpFace);
OSflux=findOSflux(G,rock,bc,interpFace);
%% nonlinear TPFA
sntp=FlowNTPFA(G,bc,fluid,Wtp,OSflux,15e6*ones(G.cells.num,1),1e-7,1000);
figure,plotCellData(G,sntp.pressure/1e6);plotWell(G,Wtp);view(3)
title('Pressure-NTPFA');h=colorbar;h.Label.String='Pressure[MPa]';clear h;
Qinj=sum(sntp.wellSol(1).flux);tt=pv/Qinj;nstep=100;dt=tt/nstep;
sntp=tracerTransport_implicit(G,rock,Wtp,sntp,dt,nstep);
figure, plotCellData(G,sntp.cres(:,end));plotWell(G,Wtp);
view(3);colorbar;title('Concentration-NTPFA');
drawnow

% NTPFA from new framework doesn't converge
% mrstVerbose off
% mrstDebug on
% fluidad = initSimpleADIFluid('phases', 'W', 'mu' , mu_value, 'rho', rho_value);
% model = GenericBlackOilModel(G, rock, fluidad, 'water', true, 'oil', false, 'gas', false);
% model = model.validateModel();
% model.FluxDiscretization.PermeabilityPotentialGradient.PermeabilityGradientDiscretization = NonlinearTwoPointFluxApproximation(model);
% model.FacilityModel = model.FacilityModel.setupWells(Wtp);
% schedule = simpleSchedule(1.0, 'W', Wtp);
% [wellSols, states] = simulateScheduleAD(state0, model, schedule);
% state = states{1};
% sntpmrst = state;
% sntpmrst = tracerTransport_implicit(G,rock,Wtp,sntpmrst,dt,nstep);


%% nonlinear MPFA
disp('nonlinear MPFA')
%snmp=NewtonMPFA(G,fluid,Wtp,OSflux,1e-9,100);
snmp=FlowNMPFA(G,bc,fluid,Wtp,OSflux,1e-9,100);
figure,plotCellData(G,snmp.pressure/1e6);plotWell(G,Wtp);axis image;view(3)
title('Pressure-NMPFA');h=colorbar;h.Label.String='Pressure[MPa]';clear h;
Qinj=sum(snmp.wellSol(1).flux);tt=1.5*pv/Qinj;nstep=100;dt=tt/nstep;
snmp=tracerTransport_implicit(G,rock,Wtp,snmp,dt,nstep);
figure, plotCellData(G,snmp.cres(:,end));plotWell(G,Wtp);
view(3);axis image;colorbar;title('Concentration-NMPFA');

% MPFA Xavier style (memory)
%disp('mpfa xavier')
%mpfastruct = computeNeumannMultiPointTrans(G, rock);
%smpfa = incompMPFA3(G, mpfastruct, Wtp);
%smpfa = tracerTransport_implicit(G,rock,Wtp,smpfa,dt,nstep);

%%
data1=abs(stp.pressure-sm.pressure)./1e6;
data2=abs(sntp.pressure-sm.pressure)./1e6;
xmin=min([data1;data2]);xmax=max([data1;data2]);
figure,subplot(1,2,1);
plotCellData(G,data1,'edgealpha',0.2);view(-80,50);axis tight off;
h=colorbar('horiz');h.Label.String='Pressure[MPa]';h.Label.FontSize=15;
title('\mid\itp_{\rmtpfa}-\itp_{\rmmfd}\rm\mid','fontsize',20);caxis([xmin xmax]);
subplot(1,2,2),plotCellData(G,data2,'edgealpha',0.2);view(-80,50);axis tight off;
h=colorbar('horiz');h.Label.String='Pressure[MPa]';h.Label.FontSize=15;
title('\mid\it{p}_{\rmntpfa}-\it{p}_{\rmmfd}\rm\mid','fontsize',20);caxis([xmin xmax]);
clear data1 data2 xmin xmax h

data=abs(stp.cres(:,end)-sm.cres(:,end));
figure,subplot(1,2,1),
plotCellData(G,data,'edgealpha',0.2);view(-80,50);axis tight off
h=colorbar('horiz');h.Label.String='Concentration';h.Label.FontSize=15;caxis([0 1]);
title('\mid\itC_{\rmtpfa}-\itC_{\rmmfd}\rm\mid','fontsize',20);
data=abs(sntp.cres(:,end)-sm.cres(:,end));
subplot(1,2,2),plotCellData(G,data,'edgealpha',0.2);view(-80,50);axis tight off
h=colorbar('horiz');h.Label.String='Concentration';h.Label.FontSize=15;caxis([0 1]);
title('\mid\itC_{\rmntpfa}-\itC_{\rmmfd}\rm\mid','fontsize',20);clear data h

figure, plot((1:nstep)'/nstep,...
             [stp.cwell ...
              stpmrst.cwell ...
              sm.cwell ...
              sntp.cwell ...
              snmp.cwell],'linewidth',2);
legend({'TPFA',...
        'mrst TPFA', ...
        'mrst MFD',...
        'NTPFA', ...
        'NMPFA'});
xlabel('pore volume injected');
ylabel('concentration');

