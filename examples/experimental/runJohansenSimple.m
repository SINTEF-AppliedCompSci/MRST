function [Gt,sol, sreport]=runJohansenSimple(method,case_name,varargin)
%% Vertical-Averaged Simulation: SLEIPNER
% Sleipner is a comercial CO2 storage site in the North Sea, where CO2 has
% been injected since 1996.
%
% In this example we simulate injection and migration of CO2 on sleipner
% using the data provided in the paper:
%
%     "Reservoir Modeling of CO2 Plume Behavior Calibrated Agains
%     Monitoring Data From Sleipner, Norway", SPE 134891
%
% The data set is avaliable online on:
%   http://www.ieaghg.org/index.php?/2009112025/modelling-network.html
%
% Provided with the paper is injection rates for 11 years. We inject 30
% years using the last injection rate the last 19 years.
%
% We demonstrate the use of C/C++-accelerated MATLAB, using the functions
%
% * processgrid (replaces processGRDECL)
% * mcomputegeometry (replaces computeGeometry)
% * mtransportVE (replaces explicitTransportVE)
%
% The last mentioned function requires that you have built the solver in
% the src/VEmex directory.

%% Display header
%clc
opt = struct('Verbose', mrstVerbose,'plot',true);
opt = merge_options(opt, varargin{:});
disp('================================================================');
disp('   Vertical averaging applied to the Sleipner model');
disp('   using C++ accelleration in the transport solver');
disp('================================================================');
disp(' ');
%a=figure();
a=40;
if(a>20)
    disp('Figure handle large close all')
   close all;
   figure() 
else
   close(a); 
end
%% Construct stratigraphic, petrophysical, and VE models
% The 3D model consists of a grid (G) and petrophysical parameters (rock).
% The VE model consists of a top-surface grid (Gt), petrophysical data
% (rock2D), and indices to the boundarcy cells where we will supply
% pressure boundary conditions. Called with a true flag, the routine will
% use C-accelerated MATLAB routines to process the data input and compute
% geometry. Once the models are created, they are stored in a data file for
% faster access at a later time.
%[Gt, rock, rock2D, bcIxVE] = makeSleipnerVEmodelSimple(true);
[G, Gt, rock, rock2D, bcIxVE] = makeJohansenVEgrid();

Gt.grav_pressure = @(g, omega) gravPressureVE_s(g, omega);
%Gt.primitives = @primitivesMimeticVE_s;
%% solver
% Select transport solver
% Use C++ acceleration if it exists - NB: requires the VEmex module
% Notice that the two solvers determine the time steps differently and
% may therefore give slightly different answers.
try
   mtransportVE;
   cpp_accel = true;
catch me%#ok
   disp('mex-file for C++ acceleration not found');
   disp(['See ', fullfile(VEROOTDIR,'src/VEmex','README'), ' for building instructions']);
   disp('Using matlab VE-transport');
   cpp_accel = false;
end
% there is isures with pressure
transport_methods={'explicit_incomp_mim','implicit_incomp_tpf','explicit_incomp_tpf',...
                   'adi_simple',...
                   'adi_OGD_simple_inst',...
                   'adi_OGD_simple_mix',...                   
                   };

%method=5;
transport_method=transport_methods{method};
[transport_solver fluidVE_h fluidVE_s, fluidADI]=...
    makeTransportSolver(transport_method,Gt,rock, rock2D,cpp_accel);



%% Set time and fluid parameters
% Fluid data are taken from paper SPE 134891
gravity on
T          = 500*year();
stopInject = 100*year();
dT         = 0.1*year();
dT_max     = 4*year();
dTpost_min = 0.1*year();
dTpost_max = 20*year();
%dTplot=1*year();
dTplot=dT_max;



%% Set well and boundary conditions
% The well is placed near the actual injection site. The injection rates
% are taken from the beforementioned paper. Hydrostatic boundary conditions
% are specified on all outer boundaries.
disp(' -> Setting well and boundary conditions');

% Set well in 3D model
wellIx = [51, 51, 6, 6];
rates = 1.4e4*meter^3/day;
%W = verticalWell([], Gt.parent, rock, wellIx(1), wellIx(2), wellIx(3):wellIx(4),...
%   'Type', 'rate', 'Val', rate, 'Radius', 0.1, 'comp_i', [1,0], 'name', 'I');


if(method<0)
  W      = verticalWell([], Gt.parent, rock, wellIx(1), wellIx(2), ...
                      1, 'Type', 'rate', 'Val', rates(1), ...
                      'Radius', 0.1, 'comp_i', [1,0], 'name', 'I');  
    
else
W      = verticalWell([], Gt.parent, rock, wellIx(1), wellIx(2), ...
                      1:Gt.parent.cartDims(3), 'Type', 'rate', 'Val', rates(1), ...
                      'Radius', 0.1, 'comp_i', [1,0], 'name', 'I','InnerProduct','ip_tpf');
end

% Well in 2D model
WVE = convertwellsVE(W, Gt.parent, Gt, rock2D);

% BC in 2D model
bcVE = addBC([], bcIxVE, 'pressure', ...
            Gt.faces.z(bcIxVE)*fluidVE_h.rho(2)*norm(gravity));
bcVE.sat = zeros(size(bcVE.face));
bcVE.h = zeros(size(bcVE.face));


%% Prepare simulations
% Compute inner products and instantiate solution structure
%{
disp(' -> Initialising solvers');
SVE = computeMimeticIPVE(Gt, rock2D, 'Innerproduct','ip_simple');
preComp = initTransportVE(Gt, rock2D);
%}
sol = initResSolVE(Gt, 300*barsa(), 0);
sol.pressure=fluidVE_h.rho(2)*norm(gravity()).*Gt.cells.z;
sol.wellSol = initWellSol(W, 300*barsa());
sol.s = height2SatVE(sol, Gt, fluidVE_h);
sol.rs=zeros(Gt.cells.num,1);
% Find trapping structure in grid. Used for calculation of trapped volumes
ts=findTrappingStructure(Gt);


%% Prepare plotting
% We will make a composite plot that consists of several parts: a 3D plot
% of the plume, a pie chart of trapped versus free volume, a plane view of
% the plume from above, and two cross-sections in the x/y directions
% through the well
opts = {'slice', wellIx, 'Saxis', [0 1-fluidVE_h.sw], 'maxH', 100, ...
   'Wadd', 500, 'view', [-85 70], 'wireH', false, 'wireS', false};
if(opt.plot)
 plotPanelVESimple(Gt.parent, Gt, W, sol, 0.0, [0 0 0 0 0 0 1], fluidVE_h, fluidADI, opts{:});
end

%% Main loop
% Run the simulation using a sequential splitting with pressure and
% transport computed in separate steps. The transport solver is formulated
% with the height of the CO2 plume as the primary unknown and the relative
% height (or saturation) must therefore be reconstructed.
t = 0;
%totVol = 0.0;
totMas = 0.0;
fprintf(1,'\nSimulating %d years of injection', convertTo(stopInject,year));
fprintf(1,' and %d years of migration\n', convertTo(T-stopInject,year));
fprintf(1,'Time: %4d years', convertTo(t,year));
stime=tic;
w = WVE;
tic
%bcVE=[];
first_post=true;
sreport={};
dT_prev=dT;

while t<T
   % Advance solution: compute pressure and then transport
   t_loc = 0;
   while t_loc < dTplot
       sol = transport_solver(sol,bcVE,w,dT);       
       assert( max(sol.s(:,1))<1+eps && min(sol.s(:,1))>-eps );
       t = t + dT;
       t_loc = t_loc +dT;
       % Add the volume injected during last time step to the total volume
       % and compute trapped and free volumes
       if ~isempty(w)
           if(method<4)
               totMas = totMas + fluidVE_h.rho(1).*w.val*dT;
           else
               rhoCO2_well= fluidADI.rhoG;%.*fluidADI.bG(sol.state.wellSol.pressure);
               totMas = totMas + rhoCO2_well.*sol.state.wellSol.qGs*dT;
               masses = massesVEADI(Gt, sol, rock2D, fluidADI, fluidVE_h);
               co2mass= masses(1)+masses(3);
               if(abs(co2mass-totMas)> totMas*1e-3)
                   disp(['Mass of co2 not in domian is ',....
                       num2str(totMas-co2mass)/1e6 ,'tonn ', num2str(abs(totMas-co2mass)/totMas),'\%'])
               end
           end
           %totVol = totVol + w.val*dT;
       end
       %    masses = massesVEADI(Gt, sol, rock2D, fluidADI, fluidVE_h);
       %    co2mass= masses(1)+masses(3);
       %    assert(co2mass<(totMas*(1+1e-3)));
       %end
       if(method < 4)
           % incompressible based
           masses = massesVE(Gt, sol, rock2D, fluidVE_h, ts);
       else
           % ADI based
           masses = massesVEADI(Gt, sol, rock2D, fluidADI, fluidVE_h);
           co2mass= masses(1)+masses(3);
           masses = massesVEADI(Gt, sol, rock2D, fluidADI, fluidVE_h, ts);
           if(~(sum(masses)<totMas*(1+1e-3)))
              disp('some problem with mass balance') 
           end
           if(~(co2mass<(totMas*(1+1e-3))))
              disp('some prolblem with mass balance') 
           end
       end
       assert((co2mass-totMas)<1e-3*totMas);
       % Check if we are to stop injecting or change injection rate
       if t>= stopInject
           w  = []; dT_max = dTpost_max; dTplot = dTpost_max;
           if(first_post)
               dT_prev=dTpost_min;
               first_post=false;
           end               
       else
           ind = min(floor(t/year)+1, numel(rates));
           rate = rates(ind);
           w.val = rate;
       end
       
       % Plotting
       
       dT=min(dT_prev*1.4,dT_max);
       if(dT<(dTplot-t_loc))
           dT_prev=dT;
       else
           if((dTplot-t_loc)>day())
            dT=dTplot-t_loc;
           end
       end
       if(t< stopInject)
          dT=min(dT,stopInject-t);
       end
   end
   if(t< stopInject)
          dTplot=min(dTplot,stopInject-t);
   end
   
   sreport{end+1}=struct('t',t,'W',W,'bcVE',bcVE,'sol',sol,'masses',[masses totMas]);%#ok
   fprintf(1,'%4d years\n', convertTo(t,year));
   %if ~(mod(t,dTplot)< dT) && t<T,
   %   continue
   %else
   if(opt.plot)
    plotPanelVESimple(Gt.parent, Gt, W, sol, t, [masses totMas],fluidVE_h , fluidADI, opts{:});
   end
   drawnow
   %end
   dTplot=min(dTplot,T-t);
end
%for i=2:numel(sreport)
%   sreport{i}.sol=rmfield(sreport{i}.sol,'twophaseJacobian')
%end
cputime=toc(stime);%#ok
save([case_name,'_',num2str(method)],'sreport','Gt','fluidVE_h','fluidVE_s','fluidADI','opts','cputime','WVE')
fprintf(1,'\n\n');
% delete C++ simulator
if cpp_accel, mtransportVE; end
%stime = toc;
disp(['Elapsed simulation time: ', num2str(toc(stime)), ' seconds.']);

%% Show the structural trap
% After the simulation has completed, we are interested in seeing how the
% location of the CO2 plume after a long migration time corresponds to the
% trapping estimates produced by trapAnalysis. This is done by finding the
% trap index of the well injection cell and then plotting the trap along
% with the final CO2 plume.

% Generate traps and find the trap corresponding to the well cells
res = trapAnalysis(Gt, false);
trap = res.trap_regions([WVE.cells]);

% Plot the areas with any significant CO2 height along with the trap in red
figure()
clf();
plotCellData(Gt, sol.h, sol.h > 0.01)
plotGrid(Gt, res.traps == trap, 'FaceColor', 'red', 'EdgeColor', 'w')
plotGrid(Gt, 'FaceColor', 'None', 'EdgeAlpha', .1);

legend({'CO2 Plume', 'Trap'})
axis tight off
view(-20, 75)
title('End of simulation CO2 compared to algorithmically determined trap')

%displayEndOfDemoMessage(mfilename)
%%
end
function [transport_solver fluidVE_h fluidVE_s fluidADI]= makeTransportSolver(solver,Gt,rock, rock2D, cpp_accel)
    mu= [6e-2*milli 8e-4]*Pascal*second;
    rho= [760 1200] .* kilogram/meter^3;
    sr= 0.21;sw= 0.11;kwm= [0.75 0.54];
    %% differnt fluids
    % Fluid for original VE formulation
    fluidVE_h    = initVEFluidHForm(Gt, 'mu' , mu, ...
                         'rho', rho, ...
                         'sr', sr, 'sw', sw, 'kwm', kwm);
    % fluid for standard incomp MRST solvers                 
    fluidVE_s = initSimpleVEFluid_s('mu' , mu , 'rho', rho, ...
                                'height'  , Gt.cells.H,...
                                'sr', [sr, sw],'kwm',kwm);                                                        
    fluidVE_s.sr = sr;
    fluidVE_s.sw = sw;
    % fluid for ad solvers gas represent CO2 oil is water
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
    if(true)
        fluidADI.bO = @(p) 1+(4.3e-5/barsa)*(p-100*barsa);
        fluidADI.BO = @(p) 1./fluidADI.bO(p);
        %p_ref=mean(Gt.cells.z).*norm(gravity).*rho(2);%average hydrostatic pressure
        T_ref=mean(Gt.cells.z) *30/1e3+273+4;% average temprature
        %[fluidADI.bG fluidADI.rhoCO2] =  boCO2(p_ref, T_ref);
        fluidADI.bG  =  boCO2(T_ref, fluidADI.rhoG);
        fluidADI.BG = @(p) 1./fluidADI.bG(p);
    end
    switch solver
        case 'explicit_incomp_mim'            
            disp(' -> Initialising solvers');
            SVE = computeMimeticIPVE(Gt, rock2D, 'Innerproduct','ip_tpf');
            preComp = initTransportVE(Gt, rock2D);
            transport_solver = @(sol,bcVE,w,dT)...
                  transport_solve_ex_mim(dT, sol, Gt, SVE,...
                                  rock, fluidVE_h, bcVE, w, preComp,...
                                  cpp_accel);
         case 'explicit_incomp_tpf'
            T = computeTrans(Gt, rock2D);
            T = T.*Gt.cells.H(gridCellNo(Gt));      
           preComp = initTransportVE(Gt, rock2D);
           transport_solver = @(sol,bcVE,w,dT)...
               transport_solve_ex_tpf(dT, sol, Gt, T,...
                                 rock2D, fluidVE_h, bcVE, w, preComp,...
                                 cpp_accel, fluidVE_s);                    
        case 'implicit_incomp_tpf'                        
            T = computeTrans(Gt, rock2D);
            T = T.*Gt.cells.H(gridCellNo(Gt));                        
            transport_solver = @(sol,bcVE_s,WVE_s,dT)...
                transport_solve_imp_tpf(dT, sol, Gt, T, fluidVE_s, WVE_s, bcVE_s, rock2D);
        case 'explicit_incomp_tpf_s'                        
            T = computeTrans(Gt, rock2D);
            T = T.*Gt.cells.H(gridCellNo(Gt));                        
            transport_solver = @(sol,bcVE_s,WVE_s,dT)...
                transport_solve_imp_tpf_s(dT, sol, Gt, T, fluidVE_s, WVE_s, bcVE_s, rock2D);           
            
        case 'adi_simple'
           s=setupSimCompVe(Gt,rock2D);
           fluidADI = addVERelperm(fluidADI, Gt,'res_oil',sw,'res_gas',sr);
           systemOG = initADISystemVE({'Oil', 'Gas'}, Gt, rock2D, fluidADI,...
               'simComponents',s,'VE',true);

           transport_solver = @(sol,bcVE_s,WVE_s,dT)...
               transport_adiOG_simple(sol, Gt, systemOG, bcVE_s, WVE_s, dT,fluidADI);
        case {'adi_OGD_simple_inst','adi_OGD_simple_mix'}           
            pressure_case='simple_instant';
            %pressure_case='simple_time_mix';
            dis_max=0.02;
            fluidADI.dis_max = dis_max;
            switch solver
                case 'adi_OGD_simple_mix'
                    %l_fac=1e-10;
                    %fluidADI.bO=@(po,rs,flag,varargin) (po-200*barsa)*l_fac+1;
                    %fluidADI.BO=@(po,rs,flag,varargin) 1./fluid.bO(po,rs,flag,varargin);
                    dis_par=0.1;% meter per year;
                    fluidADI.dis_rate=dis_par*dis_max/year;
                    %fluid.dis_rate=5e-13*H;
                    fluidADI.bO=@(po,rs,flag,varargin) fluidADI.bO(po);%(po-200*barsa)*l_fac+1;
                    fluidADI.BO=@(po,rs,flag,varargin) 1./fluidADI.bO(po,rs,flag,varargin);
                case 'adi_OGD_simple_inst'
                    %l_fac=0;
                    fluidADI.bO=@(po,rs,flag,varargin) fluidADI.bO(po);%(po-200*barsa)*l_fac+1;
                    fluidADI.BO=@(po,rs,flag,varargin) 1./fluidADI.bO(po,rs,flag,varargin);
                otherwise
            end
            fluidADI.muO=@(po,rs,flag,varargin) fluidADI.muO(po);
            fluidADI.rsSat=@(po,rs,flag,varargin)   (po*0+1)*dis_max;
            s=setupSimCompVe(Gt,rock2D);
            % important that add relperm is added after fluid properiece
            % since it is bounded to the  density
            fluidADI = addVERelperm(fluidADI, Gt,'res_oil',sw,'res_gas',sr);
            systemOG = initADISystemVE({'Oil', 'Gas','DisGas'}, Gt, rock2D, fluidADI,...
                'simComponents',s,'VE',true,'tol',1e-5);
                
            if(strcmp(pressure_case,'simple_instant'))
                systemOG.nonlinear.maxIterations=10;
            else
                systemOG.nonlinear.maxIterations=10;
            end                    
           
           transport_solver = @(sol,bcVE_s,WVE_s,dT)...
               transport_adiOGD_simple(sol, Gt, systemOG, bcVE_s, WVE_s, dT,fluidADI);     
           
           
           
         otherwise
            error('No such solver')
           transport_solver=[]; %#ok
            
    end

    
end
function sol = transport_adiOGD_simple(sol, Gt, systemOG, bcVE_s, WVE_s, dT,fluidADI)     
    
     %state=sol;
     %rmfield(sol,'state')
     %state=rmfield(state,'wellSol');
     if(isfield(sol,'state'))
         state=sol.state;
     else
        state.pressure = sol.pressure;   
        state.smax=[1-sol.extSat(:,1),sol.extSat(:,2)];
        state.smin=[1-sol.extSat(:,2),sol.extSat(:,1)];
        state.s=[1-sol.s,sol.s];
        state.rs=sol.rs;
        if(isfield(sol,'sGmax'))
            state.sGmax= sol.sGmax;
        else
            state.sGmax= state.smax(:,2);
        end 
     end
     if(isempty(WVE_s))
         W=[];
         %W = addWell([], Gt, struct('perm',ones(Gt.cells.num,1)) , [], 'Val', 0, 'Type', 'rate', 'sign', 1);
         %W.bhpLimit = 0;         
     else
         W=WVE_s;
         for i=1:numel(W)
            W(i).compi=[0,W(i).compi([2,1])]; 
         end
     end
    [state, its] = solvefiADI(state, dT, W, Gt, systemOG,'bc',bcVE_s);%#ok
    %[state, its] = solvefiADI(state, dT, W, Gt, systemOG,'bc',[]);%#ok
    sat=state.s(:,2);
    sol.s=sat;
    sol.extSat= [min(sat,sol.extSat(:,1)),max(sat,sol.extSat(:,2))];
    sol.pressure=state.pressure;
    %[s h hm] = normalizeValuesVE(Gt, sol, fluidVE_s);%#ok
    
    smax=state.sGmax;
    p=state.pressure;
    pc=fluidADI.pcOG(state.s(:,2), p,'sGmax',smax);
    pcmax=fluidADI.pcOG(smax, p,'sGmax',smax);
    drho=norm(gravity)*(fluidADI.rhoOS.*fluidADI.bO(p)-fluidADI.rhoGS.*fluidADI.bG(p));   
    h=pc./drho;
    h_max=pcmax./drho;
    h_tol=1e-1;
    if(any(h>=Gt.cells.H+h_tol) || any(h_max>=Gt.cells.H+h_tol))
       disp('Inconsistent height') 
    end
    %assert(all(h<=Gt.cells.H+1e-2));
    %assert(all(h_max<=Gt.cells.H+1e-2));
    sol.h = h;
    % due to disolution the saturation can decrease which means sGMax in
    % calculateion based on pc is wrong.
    h_max(sol.h<=0)=state.s(sol.h<=0,2).*Gt.cells.H(sol.h<=0)/(fluidADI.res_gas);
    sol.h_max = h_max;
    %sol.h_max(sol.h<=0)=state.sGmax(sol.h<=0)/(fluidADI.res_gas);
    
    %sol.h_max(sol.h<=0)=state.s(sol.h<=0,2).*Gt.cells.H(sol.h<=0)/(fluidADI.res_gas);
    s_tmp=(h.*(1-fluidADI.res_oil)+(h_max-h).*fluidADI.res_gas)./Gt.cells.H;
    assert(all(abs(s_tmp-state.s(:,2))<1e-3))
    sol.rs=state.rs;
    %Gt.cells.H.*(1-sol.s).*sol.rs/fluidADI.dis_max;
    rs_diff=(Gt.cells.H.*(1-sol.s).*sol.rs)...% all rs
        -(1-fluidADI.res_gas).*(h_max-h).*fluidADI.dis_max...% rs in the oil/water sone
        -(fluidADI.res_oil).*h.*fluidADI.dis_max;...%
    sol.rsH=rs_diff/fluidADI.dis_max+h_max;
    sol.state = state;
    sol.sGmax=state.sGmax;
end

function sol = transport_adiOG_simple(sol, Gt, systemOG, bcVE_s, WVE_s, dT,fluidADI)
     %state=sol;
     state.pressure = sol.pressure;   
     state.smax=[1-sol.extSat(:,1),sol.extSat(:,2)];
     state.smin=[1-sol.extSat(:,2),sol.extSat(:,1)];
     state.s=[1-sol.s,sol.s];
     state.rs=sol.rs;
     %state=rmfield(state,'wellSol');

     if(isempty(WVE_s))
         %W = addWell([], G, rock, 1, 'Val', 0, 'Type', 'rate', 'sign', 1);
         %W.bhpLimit = 0;
         W=[];
     else
         W=WVE_s;
         for i=1:numel(W)
            W(i).compi=[0,W(i).compi([2,1])]; 
         end
     end
    [state, its] = solvefiADI(state, dT, W, Gt, systemOG,'bc',bcVE_s);%#ok
    %[state, its] = solvefiADI(state, dT, W, Gt, systemOG,'bc',[]);%#ok
    sat=state.s(:,2);
    sol.s=sat;    
    sol.extSat= [min(sat,state.smin(:,2)),max(sat,state.smax(:,2))];
    %sol.sGmax=sol.extSat(:,2);
    sol.pressure=state.pressure;
    %%
    smax=state.smax(:,2);
    p=state.pressure;
    pc=fluidADI.pcOG(state.s(:,2), p,'sGmax',smax);
    pcmax=fluidADI.pcOG(smax, p,'sGmax',smax);
    drho=norm(gravity)*(fluidADI.rhoOS.*fluidADI.bO(p)-fluidADI.rhoGS.*fluidADI.bG(p));
    h=pc./drho;
    h_max=pcmax./drho;
    h_max(h<=0)=state.s(h<=0,2).*Gt.cells.H(h<=0)/(fluidADI.res_gas);
    s_tmp=(h.*(1-fluidADI.res_oil)+(h_max-h).*fluidADI.res_gas)./Gt.cells.H;
    assert(all(abs(s_tmp-state.s(:,2))<1e-3))
    
    %j=find(abs(s_tmp-state.s(:,2))>1e-3),
    % minum s if no compressible effects
    %state.smax(j,2)*fluidADI.res_gas/(1-fluidADI.res_oil)
    %[s h hm] = normalizeValuesVE(Gt, sol, fluidVE_s);%#ok
    sol.h = h;
    sol.h_max = h_max;
    sol.state = state;
    
end

function sol = transport_solve_ex_mim(dT, sol, Gt, SVE, rock, fluidVE, bcVE, w, preComp, cpp_accel)
    sol = solveIncompFlowVE(sol, Gt, SVE, rock, fluidVE, ...
      'bc', bcVE, 'wells', w);

   if cpp_accel
      [sol.h, sol.h_max] = mtransportVE(sol, Gt, dT, rock, ...
                                          fluidVE, 'bc', bcVE, 'wells', w, ...
                                          'gravity', norm(gravity));
   else
      sol = explicitTransportVE(sol, Gt, dT, rock, fluidVE, ...
                               'bc', bcVE, 'wells', w, 'preComp', preComp, ...
                               'intVert', false);
   end
   sol.s = height2SatVE(sol, Gt, fluidVE);
end
function sol = transport_solve_ex_tpf(dT, sol, Gt, T, rock, fluidVE_h,...
                                        bcVE, w, preComp, cpp_accel,fluidVE_s)
   sol = incompTPFA(sol, Gt, T, fluidVE_s, 'wells', w, 'bc', bcVE,'pc_form','nonwetting');
   if cpp_accel
      [sol.h, sol.h_max] = mtransportVE(sol, Gt, dT, rock, ...
                                          fluidVE_h, 'bc', bcVE, 'wells', w, ...
                                          'gravity', norm(gravity));
   else
      sol = explicitTransportVE(sol, Gt, dT, rock, fluidVE_h, ...
                               'bc', bcVE, 'wells', w, 'preComp', preComp, ...
                               'intVert', false);
   end
   sat = height2SatVE(sol, Gt, fluidVE_h);
   sol.s =sat;
   sol.extSat= [min(sat,sol.extSat(:,1)),max(sat,sol.extSat(:,2))];
end
function sol = transport_solve_imp_tpf(dT, sol, Gt, T, fluid, W2D, bc, rock2D)
    sol = incompTPFA(sol, Gt, T, fluid, 'wells', W2D, 'bc', bc,'pc_form','nonwetting');
    sol = implicitTransport(sol, Gt, dT, rock2D, fluid, ...
                            'wells', W2D, 'bc', bc, 'Verbose', false,'Trans ',T);
    [s h hm] = normalizeValuesVE(Gt, sol, fluid);%#ok
    sol.h = h;
    sol.h_max = hm;                    
    %sol = 
end
function sol = explicit_solve_imp_tpf_s(dT, sol, Gt, T, fluid, W2D, bc, rock2D)
    sol = incompTPFA(sol, Gt, T, fluid, 'wells', W2D, 'bc', bc,'pc_form','nonwetting');
    sol = explicitTransport(sol, Gt, dT, rock2D, fluid, ...
                            'wells', W2D, 'bc', bc, 'Verbose', false,'Trans ',T);
    [s h hm] = normalizeValuesVE(Gt, sol, fluid);%#ok
    sol.h = h;
    sol.h_max = hm;                    
    %sol = 
end
