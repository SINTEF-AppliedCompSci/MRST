fluid=[];
l_fac=0e-6;
fluid.rhoOS=1000;
fluid.rhoGS=600;
fluid.bO=@(po,varargin) (po-200*barsa)*l_fac+1;
fluid.BO=@(po,varargin) 1./fluid.bO(po,varargin);
fluid.bG=@(pg,varargin) pg*0.0+1;
fluid.BG=@(pg,varargin) pg*0.0+1;
fluid.muO=@(po,varargin) 0.4e-3*(po*0+1);
fluid.muG=@(pg,varargin) 1e-4*(pg*0+1);

fluid_case='hystersis';
fluid_case='simple';
fluid_case='simple_cap';
%fluid_case='cap_1D_table_simple';
%fluid_case='cap_1D_table';
H=10;
sg=linspace(0,1,100)';
Gt=[];
Gt.cells=[];
Gt.cells.H=10*ones(size(sg));

switch fluid_case
    case 'simple'
       fluid.krG=@(sg,varargin) sg;
       fluid.krOG=@(so,varargin) so;
       fluid.pcOG=@(sg, p, varargin) norm(gravity)*(fluid.rhoOS.*fluid.bO(p)-fluid.rhoGS.*fluid.bG(p)).*(sg).*Gt.cells.H;
       %fluid=rmfield(fluid,'relPerm');
       res_gas=0;
    case 'simple_cap'
        res_gas = 0.0;
        fluid = addVERelpermCapLinear(fluid,'res_gas',0.1,'beta',4,'cap_scale',0.3*H*10*(fluid.rhoOS-fluid.rhoGS),'H',Gt.cells.H,'kr_pressure',false);
    case 'cap_general'
        pc3D =@(s) 10*barsa./s.^2;
        [dpS_table, dpkrS_table] = makeSdp_table('pc3D',pc3D,'maxdP',H*drho*10*(4));
        fluid = addVERelpermTable('Sdp_table',Sdp_table,'dpkrS_table',dpkrS_table,'Gt',Gt);        
    case 'hystersis'
        res_gas = 0.5;
        fluid = addVERelperm(fluid,'res_oil',0,'res_gas',res_gas,'Gt',Gt);
    case 'cap_1D_table_simple'
        %
        S_tab=linspace(0,1,10)';
        kr_tab=S_tab;
         h_tab=S_tab*max(Gt.cells.H);
        h_mean=mean(Gt.cells.H);
%        table_co2_1d=struct('S',S_tab,'kr',kr_tab,'h',h_tab,'hH',h_tab*h_mean,...
%            'SH',S_tab*h_mean,'krH',kr_tab*h_mean);
        table_co2_1d=struct('SH',S_tab.*H,'krH',kr_tab.*H,'h',h_tab);
        
        S_tab_w=S_tab;
        kr_tab_w=S_tab_w;
        table_water_1d=struct('S',1-S_tab_w,'kr',kr_tab_w,'h',[]);
        fluid=addVERelperm1DTables(fluid,'mu',[6e-5 6.9e-4],'rho',[700 1020],'sr',[0 0],...
            'height',Gt.cells.H,'table_co2',table_co2_1d,'table_water',table_water_1d); 
     case 'cap_1D_table'
         drho=400;
        C=max(H)*0.4*drho*norm(gravity);
        alpha=0.5;
        beta = 3;
        samples=100;
        table_co2_1d=makeVEtables('pc3Dinv',@(p) (C./(p+C)).^(1/alpha),'is_kscaled',false,'kr3D',@(s) s.^beta,...
            'drho',drho,...
            'Gt',Gt,'samples',samples)
        S_tab=linspace(0,1,10)';
        tables_water_1d=struct('S',1-S_tab_w,'kr',kr_tab_w,'h',[]);
        fluid=addVERelperm1DTables(fluid,'mu',[6e-5 6.9e-4],'rho',[700 1020],'sr',[0 0],...
            'height',Gt.cells.H,'table_co2',table_co2_1d,'table_water',table_water_1d); 
     case 'cap_1D_table_pressure'
         drho=400;
        C=max(H)*0.4*drho*norm(gravity);
        alpha=0.5;
        beta = 3;
        samples=100;
        table_co2_1d = makeVEtables('pc3Dinv', @(p) (C./(p+C)).^(1/alpha),...
            'is_kscaled', false,'kr3D', @(s) s.^beta,...
            'drho', drho,...
            'Gt', Gt, 'samples', samples);
        S_tab=linspace(0,1,10)';
        tables_water_1d=struct('S',1-S_tab_w,'kr',kr_tab_w,'h',[]);
        fluid=addVERelperm1DTablesPressure(fluid,'mu',[6e-5 6.9e-4],'rho',[700 1020],'sr',[0 0],...
            'height',Gt.cells.H,'table_co2',table_co2_1d,'table_water',table_water_1d); 
                 drho=400;
      case 'cap_1D_table_kscaled'        
        surface_tension=100;  
        kscale_scale=sqrt(100/0.1)*surface_tension;
        C=max(H)*0.4*drho*norm(gravity)/kscale_scale;
        alpha=0.5;
        beta = 3;
        samples=100;
        table_co2_1d = makeVEtables('pc3Dinv', @(p) (C./(p+C)).^(1/alpha),...
            'is_kscaled', true,'kr3D', @(s) s.^beta,...
            'drho', drho,...
            'Gt', Gt, 'samples', samples);
        S_tab=linspace(0,1,10)';
        tables_water_1d=struct('S',1-S_tab_w,'kr',kr_tab_w,'h',[]);
        fluid=addVERelperm1DTablesPressure(fluid,'mu',[6e-5 6.9e-4],'rho',[700 1020],'sr',[0 0],...
            'height',Gt.cells.H,'table_co2',table_co2_1d,'table_water',table_water_1d,'rock',perm); 
    otherwise
       disp('Use deck as fluid')
end

%% relperm



p=200*barsa*ones(size(sg));
figure(1),clf
n=2;
subplot(n,1,1),cla
hold on
for smax=[0:0.2:1]
    plot(sg, fluid.krG(sg,'sGmax',smax.*ones(size(sg))));
end

% capillary pressure
subplot(n,1,2),cla
hold on
for smax=[0:0.2:1]
    plot(sg, fluid.pcOG(sg, p, 'sGmax',smax.*ones(size(sg)))/barsa);
end