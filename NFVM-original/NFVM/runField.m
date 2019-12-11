% run field model
clear;close all;clc
mrstModule add incomp mimetic
%% Grid and rock
grdecl=readGRDECL('NORNE.grdecl');usys=getUnitSystem('METRIC');
grdecl=convertInputUnits(grdecl,usys);
G=processGRDECL(grdecl,'Tolerance',0.05);G=G(1);G=computeGeometry(G);
rock=grdecl2Rock(grdecl,G.cells.indexMap);
is_pos=rock.perm(:,3)>0;rock.perm(~is_pos,3)=1e-6*min(rock.perm(is_pos,3));
is_pos=rock.perm(:,2)>0;rock.perm(~is_pos,2)=1e-6*min(rock.perm(is_pos,2));
is_pos=rock.perm(:,1)>0;rock.perm(~is_pos,1)=1e-6*min(rock.perm(is_pos,1));
figure, plotGrid(G);view(3);clear grdecl usys GG is_pos
% ----------------------------------------------------------------------
x=G.cells.centroids(:,1);y=G.cells.centroids(:,2);
ind=x>4.59e5|y>7.324e6;[~,~,K]=ind2sub(G.cartDims,G.cells.indexMap);
ind=ind|K<4;cells=find(ind);G=removeCells(G,cells);clear cells K
rock.perm=rock.perm(~ind,:);rock.poro=rock.poro(~ind);rock.ntg=rock.ntg(~ind);
% -------------------------------------------------------------------------------
figure,subplot(1,2,1);
plotCellData(G,log10(rock.perm(:,1)),'edgecolor','none');view(3);
c=colorbar('southoutside');c.Label.String='log_{10}(\itk_{\rmH})';
subplot(1,2,2);
plotCellData(G,log10(rock.perm(:,3)),'edgecolor','none');view(3)
c=colorbar('southoutside');c.Label.String='log_{10}(\itk_{\rmV})';clear c;
% ----------------------------------------------------------------
x=4.562e5;y=7.3205e6;
d=bsxfun(@minus,G.cells.centroids(:,1:2),[x y]);
d=sqrt(dot(d,d,2));[~,inj]=min(d);
x=4.583e5;y=7.3234e6;
d=bsxfun(@minus,G.cells.centroids(:,1:2),[x y]);
d=sqrt(dot(d,d,2));[~,pro]=min(d);
nz=G.cartDims(3);
[I,J,~]=ind2sub(G.cartDims,G.cells.indexMap(inj));
Wtp=verticalWell([],G,rock,I,J,1:nz,'type','bhp','val',20e6,'name','I','Comp_i', 1,'refDepth',min(G.cells.centroids(:,end)));
Wm=verticalWell([],G,rock,I,J,1:nz,'type','bhp','val',20e6,'name','I','Comp_i', 1,'innerproduct','ip_quasirt','refDepth',min(G.cells.centroids(:,end)));
[I,J,~]=ind2sub(G.cartDims,G.cells.indexMap(pro));
Wtp=verticalWell(Wtp,G,rock,I,J,1:nz,'type','bhp','val',10e6,'name','P','Comp_i', 1,'refDepth',min(G.cells.centroids(:,end)));
Wm=verticalWell(Wm,G,rock,I,J,1:nz,'type','bhp','val',10e6,'name','P','Comp_i', 1,'innerproduct','ip_quasirt','refDepth',min(G.cells.centroids(:,end)));
clear x y d ind inj pro I J nz

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
fluid=initSingleFluid('mu',1,'rho',1);
%% linear TPFA
stp=FlowTPFA(G,TransTPFA(G,rock,bc),fluid,Wtp);
figure,plotCellData(G,stp.pressure/1e6);plotWell(G,Wtp);view(3)
title('Pressure-TPFA');h=colorbar;h.Label.String='Pressure[MPa]';clear h;
Qinj=sum(stp.wellSol(1).flux);tt=pv/Qinj;nstep=100;dt=tt/nstep;
stp=tracerTransport_implicit(G,rock,Wtp,stp,dt,nstep);
figure, plotCellData(G,stp.cres(:,end));plotWell(G,Wtp);
view(3);colorbar;title('Concentration-TPFA');

%% MFD
sm=incompMimetic(initState(G,Wtp,0),G,computeMimeticIP(G,rock),fluid,'wells',Wtp);
figure, plotCellData(G,sm.pressure/1e6);plotWell(G,Wtp);view(3);
title('Pressure-MFD');h=colorbar;h.Label.String='Pressure[MPa]';clear h;
Qinj=sum(sm.wellSol(1).flux);tt=pv/Qinj;nstep=100;dt=tt/nstep;
sm=tracerTransport_implicit(G,rock,Wtp,sm,dt,nstep);
figure, plotCellData(G,sm.cres(:,end));plotWell(G,Wtp);
view(3);colorbar;title('Concentration-MFD');

%%
interpFace=findHAP(G,rock,bc);
disp(interpFace.percentage);
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

%% nonlinear MPFA
% snmp=NewtonMPFA(G,fluid,Wtp,OSflux,1e-9,100);
% snmp=FlowNMPFA(G,bc,fluid,Wtp,OSflux,1e-9,100);
% figure,plotCellData(G,snmp.pressure/1e6);plotWell(G,Wtp);axis image;view(3)
% title('Pressure-NMPFA');h=colorbar;h.Label.String='Pressure[MPa]';clear h;
% Qinj=sum(snmp.wellSol(1).flux);tt=1.5*pv/Qinj;nstep=100;dt=tt/nstep;
% snmp=tracerTransport_implicit(G,rock,Wtp,snmp,dt,nstep);
% figure, plotCellData(G,snmp.cres(:,end));plotWell(G,Wtp);
% view(3);axis image;colorbar;title('Concentration-NMPFA');

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

figure, plot((1:nstep)'/nstep,[stp.cwell sm.cwell sntp.cwell],'linewidth',2);
legend({'TPFA','MFD','NTPFA'});xlabel('pore volume injected');ylabel('concentration');

