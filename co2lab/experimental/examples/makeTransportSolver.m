function [transport_solver, fluidVE_h, fluidVE_s, fluidADI] = ...
    makeTransportSolver(solver, Gt, rock, rock2D, cpp_accel,opt)
    
    error(['This function is on-hold, and will need to be reviewed in a more ' ...
           'comprehensive manner.  As for now, it should not be presumed to work.']);
% 
%
% SYNOPSIS:
%   function [transport_solver, fluidVE_h, fluidVE_s, fluidADI] = makeTransportSolver(solver, Gt, rock, rock2D, cpp_accel,opt)
%
% DESCRIPTION:
%
% PARAMETERS:
%   solver        - 
%   Gt            - 
%   rock          - 
%   rock2D        - 
%   cpp_accel     - 
%   opt           - 
%
% RETURNS:
%   transport_solver - 
%   fluidVE_h        - 
%   fluidVE_s        - 
%   fluidADI         - 
%
% EXAMPLE:
%
% SEE ALSO:
%

    
    %% Uses hard-coded values for viscosity, density, and residual saturations
    mu= [6e-2*milli 8e-4]*Pascal*second;
    rho= [760 1200] .* kilogram/meter^3;
    sr= 0.21; sw= 0.11;kwm= [0.75 0.54];

    %% constructing different fluid objects

    % Fluid for original VE formulation
    fluidVE_h = initVEFluidHForm(Gt, ...
                                 'mu'  , mu  , ...
                                 'rho' , rho , ...
                                 'sr'  , sr  , ...
                                 'sw'  , sw  , ...
                                 'kwm' , kwm);

    % fluid for standard incomp MRST solvers                 
    fluidVE_s = initSimpleVEFluid_s('mu'     , mu , ...
                                    'rho'    , rho, ...
                                    'height' , Gt.cells.H,...
                                    'sr'     , [sr, sw], ...
                                    'kwm'    , kwm);
    fluidVE_s.sr = sr;
    fluidVE_s.sw = sw;

    % fluid for ad solvers gas represent CO2 oil is water
    fluidADI = initSimpleADIFluid('mu'  , [mu(2) mu(2) mu(1)], ...
                                  'rho' , [rho(2) rho(2), rho(1)],...
                                  'n'   , [1 1 1]); 
    wfields={'krW', 'krW', 'krG', 'pcOG', 'pcOW'};
    for i=1:numel(wfields)
        if (isfield(fluidADI, wfields{i}))
        fluidADI = rmfield(fluidADI, wfields{i});
        end
    end

    %% Adding (hard-coded) compressibility values for the ADI fluid 
    fluidADI.pvMultR =@(p) 1+(1e-5/barsa)*(p-100*barsa); % rock matrix
    
    if(is_compressible_fluid(solver))
        fluidADI.bW = @(p) 1+(4.3e-5/barsa)*(p-100*barsa);
        fluidADI.BW = @(p) 1./fluidADI.bW(p);
        p_ref=mean(Gt.cells.z).*norm(gravity).*rho(2);%average hydrostatic pressure
        T_ref=mean(Gt.cells.z) *30/1e3+273+4;% average temprature
        %[fluidADI.bG fluidADI.rhoCO2] =  boCO2(p_ref, T_ref);
        fluidADI.bG  =  boCO2(T_ref, fluidADI.rhoG);
        fluidADI.BG = @(p) 1./fluidADI.bG(p);
    end
    
    %% Choosing solver
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
      
      case {'adi_simple', 'adi_simple_incomp'}

        s=setupSimCompVe(Gt,rock2D);
        fluidADI = addVERelperm(fluidADI, Gt, 'res_water',sw, 'res_gas',sr);
        % The following line will assign stepfunction 'stepOG' and equation
        % set 'eqsfiBlackOilExplicitWellsOGVE'.
        systemOG = initADISystemVE({'Oil', 'Gas'}, Gt, rock2D, fluidADI,...
                                   'simComponents',s,'VE',true,'tol',1e-6);
        transport_solver = @(sol, bcVE_s, WVE_s, dT) ...
            transport_adiOG_simple(sol, Gt, systemOG, bcVE_s, WVE_s, dT, fluidADI);
      
      case {'adi_OGD_simple_inst','adi_OGD_simple_mix'}           
        pressure_case='simple_instant';
        %pressure_case='simple_time_mix';
        dis_max=0.02;
        %dis_max=0.02;
        fluidADI.dis_max = dis_max;
        
        fluidADI.bW=@(po,rs,flag,varargin) fluidADI.bW(po);%(po-200*barsa)*l_fac+1;
        fluidADI.BW=@(po,rs,flag,varargin) 1./fluidADI.bW(po,rs,flag,varargin);
        
        if strcmp(solver, 'adi_OGD_simple_mix')
            dis_par = 0.1;
            fluidADI.dis_rate = dis_par * dis_max / year;
        end

        fluidADI.muW=@(po,rs,flag,varargin) fluidADI.muW(po);
        fluidADI.rsSat=@(po,rs,flag,varargin)   (po*0+1)*dis_max;
        s=setupSimCompVe(Gt,rock2D);
        % important that add relperm is added after fluid properties, since
        % it is bound to density
        fluidADI = addVERelperm(fluidADI   , Gt , ...
                                'res_water'  , sw , ...
                                'res_gas'  , sr , ...
                                'top_trap' , opt.top_trap);
        systemOG = initADISystemVE({'Oil', 'Gas','DisGas'}, Gt, rock2D, fluidADI,...
                                   'simComponents' , s,...
                                   'VE'            , true, ...
                                   'tol'           , 1e-5);
        systemOG.getEquations = @eqsfiBlackOilExplicitWellsOGVE_new;   
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


function res = is_compressible_fluid(solver)
    
    switch solver
      case {'explicit_incomp_mim', ...
            'implicit_incomp_tpf', ...
            'explicit_incomp_tpf', ...
            'adi_simple_incomp'}
        res = false;
      case {'adi_simple', ...
            'adi_OGD_simple_inst', ...
            'adi_OGD_simple_mix'}
        res = true;
      otherwise
          error('Unrecognized case');
    end
end

% ----------------------------------------------------------------------------
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
        state.wellSol = sol.wellSol; % @@ Added by Odd to avoid error in the eqsfi
                                     %    function.  Should it be here???
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
    pc=fluidADI.pcWG(state.s(:,2), p,'sGmax',smax);
    pcmax=fluidADI.pcWG(smax, p,'sGmax',smax);
    drho=norm(gravity)*(fluidADI.rhoWS.*fluidADI.bW(p)-fluidADI.rhoGS.*fluidADI.bG(p));   
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
    sol.h(sol.h<=0)=0;
    sol.h_max = h_max;    
    %sol.h_max(sol.h<=0)=state.sGmax(sol.h<=0)/(fluidADI.res_gas);
    
    %sol.h_max(sol.h<=0)=state.s(sol.h<=0,2).*Gt.cells.H(sol.h<=0)/(fluidADI.res_gas);
    s_tmp=(h.*(1-fluidADI.res_water)+(h_max-h).*fluidADI.res_gas)./Gt.cells.H;
    assert(all(abs(s_tmp-state.s(:,2))<1e-3))
    sol.rs=state.rs;
    %Gt.cells.H.*(1-sol.s).*sol.rs/fluidADI.dis_max;
    rs_diff=(Gt.cells.H.*(1-sol.s).*sol.rs)...% all rs
        -(1-fluidADI.res_gas).*(h_max-h).*fluidADI.dis_max...% rs in the oil/water sone
        -(fluidADI.res_water).*h.*fluidADI.dis_max;...%
    assert(all(rs_diff>-1e-1))    
    sol.rsH=max(0,rs_diff)/fluidADI.dis_max+h_max;
    sol.state = state;
    sol.sGmax=state.sGmax;
end
% ----------------------------------------------------------------------------
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
    pc=fluidADI.pcWG(state.s(:,2), p,'sGmax',smax);
    pcmax=fluidADI.pcWG(smax, p,'sGmax',smax);
    drho=norm(gravity)*(fluidADI.rhoWS.*fluidADI.bW(p)-fluidADI.rhoGS.*fluidADI.bG(p));
    h=pc./drho;
    h_max=pcmax./drho;
    h_max(h<=0)=state.s(h<=0,2).*Gt.cells.H(h<=0)/(fluidADI.res_gas);
    s_tmp=(h.*(1-fluidADI.res_water)+(h_max-h).*fluidADI.res_gas)./Gt.cells.H;
    assert(all(abs(s_tmp-state.s(:,2))<1e-3))
    
    %j=find(abs(s_tmp-state.s(:,2))>1e-3),
    % minum s if no compressible effects
    %state.smax(j,2)*fluidADI.res_gas/(1-fluidADI.res_water)
    %[s h hm] = normalizeValuesVE(Gt, sol, fluidVE_s);%#ok
    sol.h = h;
    sol.h_max = h_max;
    sol.state = state;
    
end
% ----------------------------------------------------------------------------
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
   sol.s = height2Sat(sol, Gt, fluidVE);
end
% ----------------------------------------------------------------------------
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
   sat = height2Sat(sol, Gt, fluidVE_h);
   sol.s = sat;
   sol.extSat= [min(sat,sol.extSat(:,1)),max(sat,sol.extSat(:,2))];
end
% ----------------------------------------------------------------------------
function sol = transport_solve_imp_tpf(dT, sol, Gt, T, fluid, W2D, bc, rock2D)
    sol = incompTPFA(sol, Gt, T, fluid, 'wells', W2D, 'bc', bc,'pc_form','nonwetting');
    sol = implicitTransport(sol, Gt, dT, rock2D, fluid, ...
                            'wells', W2D, 'bc', bc, 'Verbose', false,'Trans',T);
    [s, h, hm] = normalizeValuesVE(Gt, sol, fluid);%#ok
    sol.h = h;
    sol.h_max = hm;                    
    %sol = 
end
% ----------------------------------------------------------------------------
function sol = explicit_solve_imp_tpf_s(dT, sol, Gt, T, fluid, W2D, bc, rock2D)
    sol = incompTPFA(sol, Gt, T, fluid, 'wells', W2D, 'bc', bc,'pc_form','nonwetting');
    sol = explicitTransport(sol, Gt, dT, rock2D, fluid, ...
                            'wells', W2D, 'bc', bc, 'Verbose', false,'Trans ',T);
    [s, h, hm] = normalizeValuesVE(Gt, sol, fluid);%#ok
    sol.h = h;
    sol.h_max = hm;                    
    %sol = 
end
