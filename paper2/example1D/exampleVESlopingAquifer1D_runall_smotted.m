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

mrstVerbose true
gravity on
if(true)
    n_fac=10;
    T_inj=50*year;dt_inj=10*year/5;
    T_mig=2000*year;dt_mig=100*year/5;
else
    n_fac=1;
    T_inj=7*year;dt_inj=1*year;
    T_mig=7*year;dt_mig=2*year;
end
%for smooth=[true]%false
surf_topos={'smooth','square','inf_rough','sinus'}
for k=1:numel(surf_topos)
    surf_topo=surf_topos{k};
for res_fluid=[false,true]
for depth=[2300,1300]
[nx,ny,nz] = deal(100*n_fac, 1, 1);     % Cells in Cartsian grid
[Lx,Ly,H]  = deal(30e3,10e3,50); % Physical dimensions of reservoir
total_time = 5*year;             % Total simulation time
nsteps     = 40;                 % Number of time steps in simulation
dt         = total_time/nsteps;  % Time step length
perm       = 1000;                % Permeability in milli darcies
phi        = 0.03;                % Porosity
%depth      = 1300+1000;               % Initial depth
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
W = createADIWell(W, G, rock,  floor(0.1*nx),     ...
                     'Type', 'rate', 'Val', 1*1e6/year, ...
                     'Radius', 0.125, 'Name', 'P1','Comp_i',[0 0 1]);
W = createADIWell(W, G, rock,  G.cartDims(1),     ...
                     'Type', 'bhp', 'Val', 300*barsa, ...
                     'Radius', 0.125, 'Name', 'P1','Sign',-1, 'Comp_i',[0 1 0]);

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
if(res_fluid)
    sr= 0.21;sw= 0.11;%kwm= [0.75 0.54];
else
    sr=0;sw=0;kwm=[1 1];
end
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
W(2).val=fluidADI.rhoOS*Gt.cells.z(W(2).cells)*norm(gravity);
% defnine relperm
fluid={}
fluidADI.surface_tension = 30e-3;
%fluid_names = makeVEFluidsForTest();
%fluid_names={fluid_names{[2,3,4,5,6]}}
ff_names={'sharp interface','linear cap.','P-scaled table','P-K-scaled table','S table'};
%ff_names_s={'sharp interface','linear cap.','P-scaled table','P-K-scaled table','S table'};
leg={};
%fluid={};
rock2D  = averageRock(rock, Gt);
%{
for i=1:numel(fluid_names)
    fluid{end+1} =  makeVEFluidsForTest(fluidADI, fluid_names{i},...
        'res_gas', sr,...
        'res_oil', sw,...
        'Gt', Gt, 'rock', rock2D);
    fluid{i}.name=ff_names{i};%fluid_names{i}; 
end
%}
opt=struct('res_gas', sr,...
        'res_oil', sw,...
        'Gt', Gt, 'rock', rock2D);
    results={};
ff_names={ff_names{1}}
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


for i=1:numel(ff_names);    
fluid=fluidADI;
fluid_case=ff_names{i}
 drho=(fluid.rhoOS-fluid.rhoGS);
switch fluid_case
    case 'simple'
       fluid.krG=@(sg,varargin) sg;
       fluid.krOG=@(so,varargin) so;
       fluid.pcOG=@(sg, p, varargin) norm(gravity)*(fluid.rhoOS.*fluid.bO(p)-fluid.rhoGS.*fluid.bG(p)).*(sg).*opt.Gt.cells.H;       
       fluid.res_gas=0;
       fluid.res_oil=0;
       fluid.invPc3D=@(p) 1-(sign(p+eps)+1)/2;
       fluid.kr3D=@(s) s;
    case 'integrated'
        fluid = addVERelpermIntegratedFluid(fluid, 'res_oil',opt.res_oil,...
                                            'res_gas',opt.res_gas,...
                                            'Gt',opt.Gt,...
                                            'kr_pressure',true,...
                                            'Gt',opt.Gt,...
                                            'int_poro',false,...
                                            'rock',opt.rock);
    
    case 'sharp interface'    
       fluid = addVERelperm(fluid, Gt, ...
                            'res_oil',opt.res_oil,...
                            'res_gas',opt.res_gas,...
                            'top_trap',      h_trap,...
                            'surf_topo',     surf_topo);    
    case 'linear cap.'
        fluid = addVERelpermCapLinear(fluid,...
                                      'res_gas',opt.res_gas,...
                                      'res_oil',opt.res_oil,...
                                      'beta',2,...
                                      'cap_scale',0.2*max(opt.Gt.cells.H)*10*(fluid.rhoOS-fluid.rhoGS),...
                                      'H',opt.Gt.cells.H,'kr_pressure',true);
        %{                          
            fluid = addVERelpermCapLinear(fluid,'res_gas',0.1,'beta',4,'cap_scale',...
                    0.3*H*10*(fluid.rhoOS-fluid.rhoGS),...
                    'H',Gt.cells.H,'kr_pressure',false);
        %}                              
     case 'S table' 
        C=max(opt.Gt.cells.H)*0.4*drho*norm(gravity);
        alpha=0.5;
        beta = 3;
        samples=100;
        table_co2_1d=makeVEtables('invPc3D',@(p) max((C./(p+C)).^(1/alpha),opt.res_oil),...
                                  'is_kscaled',false,...
                                  'kr3D',@(s) s.^beta,...
                                   'drho',drho,...
                                   'Gt',opt.Gt,...
                                   'samples',samples);
        S_tab=linspace(0,1,10)';
        S_tab_w=S_tab;
        kr_tab_w=S_tab_w;
        table_water_1d=struct('S',1-S_tab_w,'kr',kr_tab_w,'h',[]);
        fluid=addVERelperm1DTables(fluid,...
                                   'res_oil',opt.res_oil,...
                                   'res_gas',opt.res_gas,...
                                   'height',opt.Gt.cells.H,...
                                   'table_co2',table_co2_1d,...
                                   'table_water',table_water_1d);
     case 'P-scaled table'
        C=max(opt.Gt.cells.H)*0.4*drho*norm(gravity);
        alpha=0.5;
        beta = 3;
        samples=100;
        table_co2_1d = makeVEtables('invPc3D', @(p) max((C./(p+C)).^(1/alpha),opt.res_oil),...
            'is_kscaled', false,'kr3D', @(s) s.^beta,...
            'drho', drho,...
            'Gt', opt.Gt, 'samples', samples);
        S_tab=linspace(0,1,10)';
        S_tab_w=S_tab;
        kr_tab_w=S_tab_w;
        table_water_1d=struct('S',1-S_tab_w,'kr',kr_tab_w,'h',[]);
        fluid=addVERelperm1DTablesPressure(fluid,...
                                           'res_oil', opt.res_oil,...
                                           'res_gas', opt.res_gas,...
                                           'height',opt.Gt.cells.H,...
                                           'table_co2',table_co2_1d,...
                                           'table_water',table_water_1d,...
                                           'kr_pressure',true); 
      case 'P-K-scaled table'        
        %surface_tension=100;          
        kscale=sqrt(rock2D.poro./(rock2D.perm))*fluid.surface_tension;
        %kscale=sqrt(0.1/(100*milli*darcy))*fluid.surface_tension;        
        C=1;
        alpha=0.5;
        beta = 3;
        samples=100;
        table_co2_1d = makeVEtables('invPc3D', @(p) max((C./(p+C)).^(1/alpha),opt.res_oil),...
                                    'is_kscaled', true,....
                                    'kr3D', @(s) s.^beta,...
                                    'drho', drho,...
                                    'Gt', opt.Gt,...
                                    'samples', samples,'kscale',kscale);
        S_tab=linspace(0,1,10)';
        S_tab_w=S_tab;
        kr_tab_w=S_tab_w;
        table_water_1d=struct('S',1-S_tab_w,'kr',kr_tab_w,'h',[]);
        fluid=addVERelperm1DTablesPressure(fluid,...
                                           'res_oil', opt.res_oil,...
                                           'res_gas', opt.res_gas,...
                                           'height', opt.Gt.cells.H,...
                                           'table_co2',table_co2_1d,...
                                           'table_water',table_water_1d,...
                                           'rock',opt.rock,...
                                           'kr_pressure',true); 
    otherwise
       error('No such fluid case')
end
%fluid=fluid{1};
%fluid.pcOG = @(sg,p,varargin) 0*sg;
s=setupSimCompVe(Gt,rock2D);
if(true)
    systemOG = initADISystemVE({'Oil', 'Gas'}, Gt, rock2D, fluid,'simComponents',s,'VE',true);
else
    fluid.dis_rate=5e-13;
    dis_max=0.01;
    fluid.dis_max=dis_max;
    fluid.muO=@(po,rs,flag,varargin) fluidADI.muO(po);
    fluid.rsSat=@(po,rs,flag,varargin)   (po*0+1)*dis_max;
    systemOG = initADISystemVE({'Oil', 'Gas','DisGas'}, Gt, rock2D, fluid,'simComponents',s,'VE',true);
    
end
%systemOG  = initADISystemVE({'Oil', 'Gas','DisGas'}, Gt, rock, fluid)%, 'tol_mb', 0);




%% Run the schedule setup in the file
% Before we can run the schedule, we make sure that we have an initial
% hydrostatic pressure distribution. Then we pick the schedule from the
% input deck and start the simulator.
%x0 = initEclipseState(G, deck, initEclipseFluid(deck));
z  = G.cells.centroids(:,3);
clear x0;
%x0.pressure = ipress*barsa +(z(:)-z(end))*norm(gravity)*fluid.rhoOS;
x0.pressure = W(2).val +(z(:)-z(W(2).cells))*norm(gravity)*fluid.rhoOS;
x0.s(:,1)=ones(G.cells.num,1);
x0.s(:,2)=zeros(G.cells.num,1);
x0.rs=ones(G.cells.num,1)*0.0;
x0.smax=x0.s;
x0.smin=x0.s;
x0.sGmax=x0.s(:,2);
%x0.s=x0.s(:,[2,1]);
%x0.s=x0.z(:,[2,1]);
%dt=ones(5,1)*0.5*year
dt=linspace(0.1*year,dt_inj,10)';
dt_shift=sum(dt);
dt=[dt;ones(floor((T_inj-dt_shift)/dt_inj),1)*dt_inj];
dt=[dt;T_inj-sum(dt)];
dt_post=linspace(0.5*year,dt_mig,5)';%'*year;
dt_shift=sum(dt_post);
dt_post=[dt_post;ones(floor((T_mig-dt_shift)/dt_mig),1)*dt_mig];%*year];
dt_post=[dt_post;T_mig-sum(dt_post)];
%%
control=struct('W',[],'step',struct('val',[],'control',[]));
W_post=W(2);
control.W={W,W_post};
control.step.val=[dt;dt_post];
%control.step.control=ones(size(dt));
control.step.control=[ones(size(dt));ones(size(dt_post))*2];

%[wellSols, states] = runScheduleADI(x0, G, rock, systemOG, deck.SCHEDULE);
%profile off
%profile on
systemOG.nonlinear.linesearch=true;
systemOG.nonlinear.maxIterations=10;
[wellSols, states] = runMrstADI(x0, Gt, systemOG, control,'force_step',false,'dt_min', 0.5*year,'report_all',false);
results{end+1}=struct('wellSols',{wellSols},'states',{states});
%profile off
%profile viewer
end
zz=Gt.cells.z;
xc=Gt.cells.centroids(:,1)/1e3;

        if(res_fluid)
            fname=['topo_model_',surf_topo,'_res_fluid','_','_depth_',num2str(depth)];
        else
            fname=['topo_model_',surf_topo,'_nores_fluid','_','_depth_',num2str(depth)];
        end

 save(['test_data_smotted/',fname],'results','z','xc','opt')   
 
end
end
end
return
%%
xx=xc(150:650);ff=exp(-((xx-xc(400))/(0.3)).^2);plot(xx,ff)
ff=ff/sum(ff);
clf,plot(xc,s,xc,filter2(ff./sum(ff),s))
figure(1),clf,hold on

hold on;
marker={'s','*','d','o','<'}
marker={'ks','b*','md','ro','<'}
marker={'k','g','m','r','b'}
use_filter=true;
for k=1:numel(results)
    states=results{k}.states;
    p0=x0.pressure;
    xc = Gt.cells.centroids(:,1)/1e3;
    zt = Gt.cells.z;
    zb = zt + Gt.cells.H;
    ns=numel(states);
    %for nn=[numel(dt),ns-10]%ns];%[numel(dt_inj),ns]
    for nn=[ns-10,ns-10]%ns];%[numel(dt_inj),ns]    
        %for nn=ns:ns
        
        state=states{nn};
        
        if(isfield(fluid,'dis_max'))
            smax=state.sGmax;
            rsH=state.s(:,1).*state.rs/fluid.dis_max;
        else
            smax=state.smax(:,2);
            rsH=state.s(:,1)*0;
        end
        diff=max(zt)-min(zt);
        %subplot(2,2,1),cla
        %title('pressure')
        %plotCellData(G,state.pressure/barsa);colorbar('horiz'), caxis([100 600])
        
        %subplot(2,2,2),cla
        %title('saturation')
        %plotCellData(G,state.s(:,2));colorbar('horiz'), caxis([0 1])
        
        % plot as VE
        %{
    subplot(3,1,2),hold on
    plot(xc,state.pressure/barsa); set(gca,'YLim',[100 300]);
    subplot(3,1,3),hold on
    plot(xc,(state.pressure-p0)/barsa); %set(gca,'YLim',[0 200]);
    subplot(3,1,1),hold on
        %}
        sG=free_sg(state.s(:,2),smax,opt);
        if(use_filter==false)
            if(isfield(fluid,'dis_max'))
                plot(xc,state.s(:,2).*Gt.cells.H,xc,smax.*Gt.cells.H,xc,smax.*Gt.cells.H+rsH); %set(gca,'YLim',[0 20]);
            else
                %plot(xc,state.s(:,2).*Gt.cells.H,['-',marker{k}],xc,smax.*Gt.cells.H,['-',marker{k}]);
                plot(xc,sG.*Gt.cells.H,['-',marker{k}],'LineWidth',2)
                plot(xc,state.s(:,2).*Gt.cells.H,['.',marker{k}],'LineWidth',2)
                %plot(xc,smax.*Gt.cells.H,['-',marker{k}],'LineWidth',2) 
            end
        else
           plot(xc,filter2(ff,sG.*Gt.cells.H),['-',marker{k}],'LineWidth',2)
           plot(xc,filter2(ff,smax.*Gt.cells.H),['.',marker{k}],'LineWidth',2) 
           plot(xc,(sG.*Gt.cells.H))%,['-',marker{k}],'LineWidth',2)
           plot(xc,(smax.*Gt.cells.H))%,['.',marker{k}],'LineWidth',2) 
        end
        %%
    end
    %%
end
set(gca,'FontSize',19)
set(gca,'Ydir','reverse')
h1 = gca;
h3 = axes('Position',get(h1,'Position'));
x=1:2;
leg_names=ff_names;
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
    set(f2(f),'Color',marker{f});
end
set(h3,'FontSize',16);
bb=legend(leg_names,'Location','SouthEast','FontSize',16);
set(gca,'FontSize',19)


return

%%

%% Plot results
%figure
states=results{1}.states;
p0=x0.pressure;
xc = Gt.cells.centroids(:,1)/1e3;
zt = Gt.cells.z;
zb = zt + Gt.cells.H;
for nn=1:numel(states)
    figure(1),clf
    state=states{nn};
    
    if(isfield(fluid,'dis_max'))
        smax=state.sGmax;
        rsH=state.s(:,1).*state.rs/fluid.dis_max;
     else        
        smax=state.smax(:,2);
        rsH=state.s(:,1)*0;
     end
    diff=max(zt)-min(zt); 
    %subplot(2,2,1),cla
    %title('pressure')
    %plotCellData(G,state.pressure/barsa);colorbar('horiz'), caxis([100 600])
    
    %subplot(2,2,2),cla
    %title('saturation')
    %plotCellData(G,state.s(:,2));colorbar('horiz'), caxis([0 1])
    
    % plot as VE
    subplot(3,1,2),cla   
    plot(xc,state.pressure/barsa); set(gca,'YLim',[100 300]);
    subplot(3,1,3),cla   
    plot(xc,(state.pressure-p0)/barsa); %set(gca,'YLim',[0 200]);
    subplot(3,1,1),cla,hold on
    if(isfield(fluid,'dis_max'))
    plot(xc,state.s(:,2).*Gt.cells.H,xc,smax.*Gt.cells.H,xc,smax.*Gt.cells.H+rsH); %set(gca,'YLim',[0 20]);
    else
       plot(xc,state.s(:,2).*Gt.cells.H,xc,smax.*Gt.cells.H);
       hh=state.s(:,2).*Gt.cells.H;
       plot(xc,filter2(ff./sum(ff),hh))
    end
    %%
    figure(2)
     clf
     
    %
    patch(xc([1 1:end end]), [zt(end)-20; zt; zt(end)-20],.7*[1 1 1]);
    patch(xc([1 1:end end]), [zb(1)+20; zb; zb(1)+20],.7*[1 1 1]);
    patch(xc([1:end end:-1:1]), ...
      [zt + Gt.cells.H.*state.s(:,2); zt(end:-1:1)], 'g')
  %
  %%{
    patch(xc([1:end end:-1:1]), ...
      [zt + Gt.cells.H.*state.s(:,2); zt(end:-1:1)+Gt.cells.H.*(smax(end:-1:1))],[1 1 0])
    patch(xc([1:end end:-1:1]), ...
      [zt + Gt.cells.H.*smax; zt(end:-1:1)+Gt.cells.H.*(smax(end:-1:1)+rsH(end:-1:1))], [0 0.1 1])
    patch(xc([1:end end:-1:1]), ...
      [zt + Gt.cells.H.*(smax+rsH); zb(end:-1:1)], [0 0 0.6])
  %}  
    line(xc([1:end]),zt + Gt.cells.H.*smax,'LineWidth',2,'Color','k')
    patch(xc([1 1:end end]), [zt(end)-20; zt; zt(end)-20],.7*[1 1 1]);
    set(gca,'YDir','reverse'), axis tight
    %%
    drawnow;
    
    pause(0.1)
end
return
%%
plot(fluid.krG(linspace(0,1,G.cells.num)',700*barsa*ones(G.cells.num,1)))
%%
s=linspace(0,1,G.cells.num)';
s_ad=initVariablesADI(s);
krG=fluid.krG(s_ad,700*barsa*ones(G.cells.num,1));
plot(s,diag(krG.jac{1}))
%%
pp=states{1}.pressure
pp=x0.pressure
%figure(99),plot(fluidADI.bG(pp)*fluid.rhoGS,Gt.cells.z);axis([0 1000 0 3000]);set(gca,'YDir','reverse')
figure(99),plot(xc,fluidADI.bG(pp)*fluid.rhoGS);set(gca,'YDir','reverse')
