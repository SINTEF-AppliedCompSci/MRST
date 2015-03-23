%% VE simulation in a standard black-oil solver
%  In this example we show how to set up a standard format black-oil
%  model that can be used to simulate a VE model. For the actual
%  simulation,  we use the fully-implicit solver in MRST from the 'ad-fi'
%  module, which is based on automatic differentiation. 

try
   require deckformat ad-fi
catch %#ok<CTCH>
   mrstModule add deckformat ad-fi
end

%% Parameters for the simulation
%close all
mrstModule add upscaling
mrstVerbose true
gravity on
force_time_step=true;
if(true)
    n_fac=10;
    T_inj=50*year;dt_inj=10*year/5;
    T_mig=2000*year;dt_mig=100*year/5;
else
    if(false)
        n_fac=1;
        T_inj=7*year;dt_inj=1*year;
        T_mig=7*year;dt_mig=2*year;
    else
        n_fac=1;
        T_inj=50*year;dt_inj=5*year;
        T_mig=100*year;dt_mig=20*year;        
    end
end
%for smooth=[true]%false
surf_topos={'smooth','square','inf_rough','sinus'}
 surf_topo=surf_topos{end};

[nx,ny,nz] = deal(100*n_fac, 1, 1);     % Cells in Cartsian grid
[Lx,Ly,H]  = deal(30e3,10e3,50); % Physical dimensions of reservoir
total_time = 5*year;             % Total simulation time
nsteps     = 40;                 % Number of time steps in simulation
dt         = total_time/nsteps;  % Time step length
perm       = 1000;                % Permeability in milli darcies
phi        = 0.03;                % Porosity
depth      = 1300+1000;               % Initial depth
%ipress     = 300;
%dp         = 70;

%% Create input deck and construct grid
% Create an input deck that can be used together with the fully-implicit
% solver from the 'ad-fi' module. Since the grid is constructed as part of
% setting up the input deck, we obtain it directly. 
G=cartGrid([nx,ny,nz],[Lx,Ly,H]);
x=G.nodes.coords(:,1);
LL=Lx*2/3;
G_org=G;
if(~(strcmp(surf_topo,'smooth')))
    G.nodes.coords(:,3)=G.nodes.coords(:,3)+depth-LL*sin(x/LL)*tan(phi);%+2*sin(2*pi*x/0.3e3);
    G_org.nodes.coords(:,3)=G_org.nodes.coords(:,3)+depth-LL*sin(x/LL)*tan(phi)+2*sin(2*pi*x/0.3e3);
    G_org=computeGeometry(G_org);
    Gt_org= topSurfaceGrid(G_org);
else
    G.nodes.coords(:,3)=G.nodes.coords(:,3)+depth-LL*sin(x/LL)*tan(phi)+2*sin(2*pi*x/0.3e3);
end
%G.nodes.coords(:,3)=G.nodes.coords(:,3)+depth+5*sin(x*10*pi/Lx)-7*sin(x*3*pi/Lx)-x*tan(phi);
G=computeGeometry(G);
%{
[deck, G] = sinusDeckAdi_GasOilDisolved([nx ny nz], [Lx Ly H], nsteps, dt, ...
                         -.1*pi/180, depth, phi, perm, ...
                         (H*phi*Lx*Ly)*0.2*day/year, ipress,dp);
%}
rock=struct('perm',100*milli*darcy*ones(G.cells.num,1),'poro',0.2*ones(G.cells.num,1));
% Alternatively, we could read deck from file and construct the grid
% deck = readEclipseDeck( ...
%    fullfile(VEROOTDIR,'data','decks','sinusDeckAdi.DATA');
% G = initEclipseGrid(deck);
%%
W=[];
W = createSampleWell(W, G, rock,  floor(0.1*nx),     ...
                     'Type', 'rate', 'Val', 1*1e6/year, ...
                     'Radius', 0.125, 'Name', 'P1','Comp_i',[0 0 1]);
W = createSampleWell(W, G, rock,  G.cartDims(1),     ...
                     'Type', 'bhp', 'Val', 300*barsa, ...
                     'Radius', 0.125, 'Name', 'P1','Sign',-1, 'Comp_i',[0 0 1]);

figure, plotGrid(G),plotWell(G,W),view([0 -1 0]), box on
                      


%% Initialize data structures
% First, we convert the input deck to SI units, which is the unit system
% used by MRST. Second, we initialize the rock parameters from the deck;
% the resulting data structure may have to be post-processed to remove
% inactive cells. Then we set up the fluid object and tell the ad-fi solver
% that that we are working with an oil-gas system.
%deck  = convertDeckUnits(deck);
%rock  = initEclipseRock(deck);
%rock  = compressRock(rock, G.cells.indexMap);
%fluid = initDeckADIFluid(deck);
% set the capillary pressure and the VE relperms explicitely
Gt = topSurfaceGrid(G);



%Temp=mean(Gt.cells.z*30/1e3+274);
mu= [6e-2*milli 8e-4]*Pascal*second;
rho= [760 1100] .* kilogram/meter^3;
%if(res_fluid)
%    sr= 0.21;sw= 0.11;%kwm= [0.75 0.54];
%else
    sr=0;sw=0;kwm=[1 1];
%end
%sr= 0.21;sw= 0.11;kwm= [1 1];
fluidADI = initSimpleADIFluid('mu',[mu(2) mu(2) mu(1)],...
                                'rho',[rho(2) rho(2), rho(1)],...
                                'n',[1 1 1]); 
wfields={'krO''krW','krG','pcOG','pcOW'};
for i=1:numel(wfields)
    if(isfield(fluidADI,wfields{i}))
        fluidADI=rmfield(fluidADI,wfields{i});
    end
end
%%
fluidADI.pvMultR =@(p) 1+(1e-5/barsa)*(p-100*barsa);
fluidADI.bO = @(p,varargin) 1+(4.3e-5/barsa)*(p-100*barsa);
fluidADI.BO = @(p, varargin) 1./fluidADI.bO(p);
%Temp=mean(Gt.cells.z*20/1e3+273+15);
Temp=Gt.cells.z*30/1e3+274+12;
%Temp=(Gt.cells.z-1000)*10/1e3+274+40;
p0=Gt.cells.z(:)*norm(gravity)*1100;
%pp0 = W(2).val +(z(:)-z(W(2).cells))*norm(gravity)*fluid.rhoOS;
fluidADI.bG  =  boCO2(Temp, fluidADI.rhoGS);fluidADI.BG = @(p) 1./fluidADI.bG(p);
figure(44),clf,plot(fluidADI.bG(p0)*fluidADI.rhoGS,Gt.cells.z);axis([0 1000 0 3000]);set(gca,'YDir','reverse')
figure(99),clf,plot(Gt.cells.centroids(:,1)/1e3,fluidADI.bG(p0)*fluidADI.rhoGS);set(gca,'YDir','reverse')
grid on

%%
%plot(Gt.cells.z,fluidADI.bG(p0)*fluid.rhoGS)
%return
%%
%{
fluidADI.bG= @(p) 1+0*p;
fluidADI.BG=@(p)  1+0*p;
fluidADI.bO= @(p,varargin) 1+0*p;
fluidADI.BO=@(p,varargin)  1+0*p;
fluidADI.rhoGS=600;fluidADI.rhoOS=1000;
%}

% defnine relperm
fluid={}
fluidADI.surface_tension = 30e-3;
%fluid_names = makeVEFluidsForTest();
%fluid_names={fluid_names{[2,3,4,5,6]}}
ff_names={'sharp interface','linear cap.','P-scaled table','P-K-scaled table','S table'};
%ff_names_s={'sharp interface','linear cap.','P-scaled table','P-K-scaled table','S table'};
leg={};

rock2D  = averageRock(rock, Gt);

opt=struct('res_gas', sr,...
        'res_water', sw,...
        'Gt', Gt, 'rock', rock2D);
    results={};
ff_names={ff_names{3}}
if(~(strcmp(surf_topo,'smooth')))
%% make smoothed parameters
xc=Gt.cells.centroids(:,1)/1e3;
 xx=xc(150:650);ff=exp(-((xx-xc(floor(nx/2)))/(0.3)).^2);plot(xx,ff)
ff=exp(-((xx-xc(400))/(0.3)).^2);ff=ff/sum(ff);plot(ff)
%%
xc=Gt.cells.centroids(:,1)/1e3;
z=Gt.cells.z;
z_org=Gt_org.cells.z;
z_new=max(z_org)*ones(size(z_org));
%assert(z_org(1)==z_new(1));
%%
for i=1:numel(z_new)-1
    z_new(i)=max(z_org(i:end));
end
z_new(end)=max(z_new(end-1),z_org(end))
plot(xc,[z,z_org,z_new])
%plot(xc,[z_org,z_new])
legend('org','new')
shirt=100
pp=numel(xc)-700;
axis([xc(pp-shirt) xc(pp) min(z_org(pp-shirt:pp)) max(z_org(pp-shirt:pp))])
plot(xc,filter2(ff,z_new-z_org))
h_trap=filter2(ff,z_new-z_org)
%%
else
    h_trap=[];
end
%%
clear surf_topos
surf_topos={'smooth','square','inf_rough','sinus'}
fluid=cell(numel(surf_topos),1);
for i=1:numel(surf_topos) 
    surf_topo=surf_topos{i}
    fluid{i} = addVERelperm(fluidADI, opt.Gt, ...
                        'res_water',opt.res_water,...
                        'res_gas',opt.res_gas,...
                        'top_trap',      h_trap,...
                        'surf_topo',     surf_topo);  
end
%{
fluid{end+1} = addVERelperm(fluidADI, opt.Gt, ...
                        'res_water',opt.res_water,...
                        'res_gas',opt.res_gas,...
                        'top_trap',      2.*ones(size(h_trap)),...
                        'surf_topo',     surf_topo);
 %}
%%
figure(33),clf,hold on
ns=100;
sGv=linspace(0,1*4/50,ns);
krG=nan(Gt.cells.num,ns);
krO=nan(Gt.cells.num,ns);
cno=floor(1000*20/30);
muO=fluid{1}.muO(300*barsa);
muG=fluid{1}.muG(300*barsa);
markers={'k','g','m','r','b'}
for k=1:numel(fluid)
    for j=1:numel(sGv)
        sG=sGv(j)*ones(Gt.cells.num,1);
        krG(:,j)=fluid{k}.krG(sG,p0);
        krO(:,j)=fluid{k}.krO(1-sG,p0);
    end
    %fw=(krO.*krG/(muO.*muG))./((krG./muG) + (krO./muO));
    plot(sGv*50,krG(cno,:),markers{k},'LineWidth',2)
    %plot(sGv*50,fw(cno,:))
end
set(gca,'FontSize',16)
%axes(h1)
us=load('data/upscaled_relperm_theta');%us=us.krCO2
plot(us.sat_mat*50,squeeze(us.krCO2(:,:,1)),'*-')
f3=plot(us.sat_mat*50,squeeze(us.krCO2(:,:,1)),'*-')%,sat_mat,kr{2}(:,:))
axis([0 4 0 0.08])
%
set(gca,'FontSize',19)
%set(gca,'Ydir','reverse')
h1 = gca;
h3 = axes('Position',get(h1,'Position'));
x=1:2;
leg_names={'smooth','square','inf rough','sinus'};
tmp=nan(numel(x),numel(leg_names));
%tmp=repmat(x,numel(leg2),1);
f2=plot(x,tmp)
axis off;
set(gca,'XtickLabel',[])
for f=1:numel(f2)
    set(f2(f),'Marker','none');
    set(f2(f),'MarkerFaceColor','none');
    set(f2(f),'LineStyle','-');
    %set(f2(f),'Color',get(a(f),'Color'));
    %set(f2(f),'Color',mycol(f,:));
    set(f2(f),'Color',markers{f});
end
set(h3,'FontSize',16);
bb=legend(leg_names,'Location','SouthEast','FontSize',16);
set(gca,'FontSize',19)
hold on
%%
% from running simple1D_VE_upscaling.m

%%
h3 = axes('Position',get(h1,'Position'));
x=1:2;
leg_names={'\theta=0','\theta=0.0162','\theta=0.03'};
tmp=nan(numel(x),numel(leg_names));
%tmp=repmat(x,numel(leg2),1);
f2=plot(x,tmp)
axis off;
set(gca,'XtickLabel',[])
for f=1:numel(f2)
    set(f2(f),'Marker',get(f3(f),'Marker'));
    set(f2(f),'MarkerFaceColor','none');
    set(f2(f),'LineStyle','-');
    set(f2(f),'Color',get(f3(f),'Color'));
    %set(f2(f),'Color',mycol(f,:));
    %set(f2(f),'Color',markers{f});
end
set(h3,'FontSize',16);
bb=legend(leg_names,'Location','NorthWest','FontSize',16);
%set(gca,'FontSize',19)
axes(h3)
set(gcf,'PaperPositionMode','auto')
print('-depsc2','figs/upscaled_relperms.eps')
%dip=diff(Gt.cells.z)/diff(Gt.cells.centroids(:,1));
%figure(1),plot(dip)
