moduleCheck('co2lab','ad-fi','ad-props','co2lab-ad');
%%
% parameters for plotting
param=[0.0,0.2;...
    0.0,0.0;...
    0.2,0.2;...
    0.2, 0.0;...
    0.21, 0.11...
    ];
for kk=1:size(param,1)
    src=param(kk,1); srw=param(kk,2);
    gravity on;
    T_ref=277+1000*30/1e3;
    mu= [6e-2*milli 8e-4]*Pascal*second;
    rho= [760 1200] .* kilogram/meter^3;
    %src= 0.21;srw= 0.11;
    kwm= [0.75 0.54];    %src= 0.0;srw= 0.2;
    plot_hyst=(src>0);
    %myname=
    myname=['resc_',num2str(src*100),'_resw_',num2str(srw*100)];
    H=25;
    sample=100;
    G=cartGrid([sample 1 1],[100 100 H]);
    G=computeGeometry(G);
    Gt=topSurfaceGrid(G);
    rock=struct('perm',1000*milli*darcy*ones(G.cells.num,1),'poro',0.4*ones(G.cells.num,1));
    rock.parent=rock;
    fluidADI = initSimpleADIFluid('mu',[mu(2) mu(2) mu(1)],...
                                  'rho',[rho(2) rho(2), rho(1)],...
                                  'n',[1 1 1]);
    wfields={'krW', 'krO','krG','pcOG','pcOW'};
    for i=1:numel(wfields)
        if(isfield(fluidADI,wfields{i}))
            fluidADI=rmfield(fluidADI,wfields{i});
        end
    end
    fluidADI.pvMultR =@(p) 1+(1e-5/barsa)*(p-100*barsa);
    fluidADI.bW = @(p) 1+(4.3e-5/barsa)*(p-100*barsa);
    fluidADI.BW = @(p) 1./fluidADI.bW(p);
    
    fluidADI.bG  =  boCO2(T_ref, fluidADI.rhoGS);
    fluidADI.BG = @(p) 1./fluidADI.bG(p);
    
    fluidADI.bG= @(p) 1+0*p;fluidADI.BG=@(p)  1+0*p;fluidADI.bW= @(p) 1+0*p;fluidADI.BW=@(p)  1+0*p;
    fluidADI.rhoGS=600;fluidADI.rhoWS=1000;
    
    
    clear fluid
    fluid{1}=addVERelperm(fluidADI, Gt, ...
        'res_water',       srw,...
        'res_gas',       src,...
        'kro',           1,...
        'krg',           1,...
        'top_trap',      [],...
        'surf_topo',     'smooth');
    fluid{1}.name='linear';
    fluid{2}=addVERelperm(fluidADI, Gt, ...
        'res_water',       srw,...
        'res_gas',       src,...
        'kro',           1,...
        'krg',           1,...
        'top_trap',      5,...
        'surf_topo',     'square');
    fluid{2}.name='square';%#ok
    %%
    fluid={};
    fluidADI.surface_tension = 30e-3;
    fluid_names = makeVEFluidsForTest();
    fluid_names={fluid_names{[2,3,4,5,6]}};%#ok
    ff_names={'sharp interface','linear cap.','P-scaled table','P-K-scaled table','S table'};
    leg={};
    %fluid={};
    for i=1:numel(fluid_names)
        fluid{end+1} =  makeVEFluidsForTest(fluidADI, fluid_names{i},...
            'res_gas', 0,...
            'res_water', 0,...
            'Gt', Gt, 'rock', rock,'res_water',srw,'res_gas',src);%#ok
        fluid{i}.name=ff_names{i};%#ok 
    end
    
    %%
    
    s=linspace(0,1,sample)';
    p=100*barsa*ones(size(s));
    sco_max=0.6*ones(size(s));
    figure(1),clf,hold on
    
    for i=1:numel(fluid)
        subplot(2,1,1),hold on
        plot(s,fluid{i}.krG(s,p,'sGmax',s),'b',s,fluid{i}.krG(s,p,'sGmax',sco_max),'r')
        subplot(2,1,2),hold on
        plot(s,fluid{i}.pcWG(s,p,'sGmax',s),'b',s,fluid{i}.pcWG(s,p,'sGmax',sco_max),'r')
    end
    %%
    krG=[];krGm=[];pcWG=[];pcWGm=[];
    for i=1:numel(fluid)
        krG=[krG, fluid{i}.krG(s,p,'sGmax',s)+i*1e-3];%#ok
        krGm=[krGm, fluid{i}.krG(s,p,'sGmax',sco_max)];%#ok
        pcWG=[pcWG,fluid{i}.pcWG(s,p,'sGmax',s)];%#ok
        pcWGm=[pcWGm,fluid{i}.pcWG(s,p,'sGmax',sco_max)];%#ok
        leg{i}=fluid{i}.name;%#ok
    end
    %%
    figure(11),clf,hold on;
    plot(s,krG,'-','LineWidth',2);
    legend(leg{:},'Location','NorthWest')
    if(plot_hyst)
        plot(s,krGm,'--','LineWidth',2)
    end
    box on
    set(gca,'FontSize',16)
    set(gca,'LineWidth',2)
    %my
    mkdir('figs/relperm')
    print(['figs/relperm/',myname,'_krG'])
    %%
    dpH=9.8*400*H/barsa;
    figure(12),clf,hold on;
    plot(s,pcWG/barsa,'-','LineWidth',2)
    legend(leg{:},'Location','NorthWest')
    if(plot_hyst)
        plot(s,pcWGm/barsa,'--','LineWidth',2)
    end
    box on
    line([0 1],[dpH dpH],'LineWidth',4,'Color','k')
    axis([0 1 0 2*dpH])
    set(gca,'FontSize',16)
    set(gca,'LineWidth',2)
    %my
    print(['figs/relperm/',myname,'_pcWG'])
    %%
    dpH=9.8*400*H;
    figure(13),clf,hold on;
    plot(s,pcWG/dpH,'-','LineWidth',2)
    legend(leg{:},'Location','NorthWest')
    if(plot_hyst)
        plot(s,pcWGm/dpH,'--','LineWidth',2)
    end
    box on
    line([0 1],[1 1],'LineWidth',4,'Color','k')
    axis([0 1 0 2])
    set(gca,'FontSize',16)
    set(gca,'LineWidth',2)
    %my
    print(['figs/relperm/',myname,'_hOG'])
    
    %%
    sco_max=mean(sco_max);
    
    %%
    figure(2),clf,hold on
    kscale=sqrt(rock.poro./(rock.perm))*fluidADI.surface_tension;
    if(plot_hyst)
        h_int=0.5*sco_max*H;
        h_int_max=sco_max*H;
    else
        h_int=sco_max*H;
        h_int_max=sco_max*H;
    end
    leg={};
    ss=[];hh=[];krk=[];ss_max=[];hh_max=[];
    for i=1:numel(fluid)
        [s,pc,kr,SH,krH, s_max,fval]= veRelpermTester(h_int, p, fluid{i}, H, 'samples',100,'hs_max',h_int_max,'kscale',kscale);%#ok
        ss=[ss,1-fval.s_h];%#ok
        ss_max=[ss_max,1-fval.s_hmax];%#ok
        hh_max=[hh_max,h_int_max-fval.h_max];%#ok
        krk=[krk,fval.kr_h];%#ok
        hh=[hh,h_int-fval.h];%#ok
        %plot(1-fval.s_h,h_int-fval.h,'LineWidth',2)
        %set(gca,'YDir','reverse')
        %axis([0 1,0 50])
        leg{i}=fluid{i}.name;%#ok
    end
    plot(ss,hh,'LineWidth',2);
    if(plot_hyst)
        plot(ss_max,hh_max,'--','LineWidth',2)
    end
    legend(leg{:},'Location','NorthWest')
    set(gca,'YDir','reverse')
    axis([0 1.02,0 H])
    box on
    set(gca,'LineWidth',2)
    set(gca,'FontSize',16)
    %my
    print(['figs/relperm/',myname,'_sz_',num2str(floor(max(h_int)))])
end
%{
 figure(3)
 plot(krk,hh,'LineWidth',2)
 legend(leg{:},'Location','SouthEast')
 set(gca,'YDir','reverse')
 axis([0 1,0 H])
 box on
 set(gca,'LineWidth',2)
%}