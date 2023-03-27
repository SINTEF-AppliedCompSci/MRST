function solver = solverFactoryVE(Gt, fluid, rock, system, varargin)
   opt = struct('transport',      'implicit',...
                'formulation',    's', ...
                'transportTolerance', 1e-3, ...
                'discretization', 'tpfa' ...
                );
   
   opt = merge_options(opt, varargin{:});
   
   fullyImplicit = ~isempty(system);
   useS = strcmpi(opt.formulation, 's');
   
   % Misc constants and common precomp
   rock2D    = averageRock(rock, Gt);
   T = [];
   
   if fullyImplicit
       require ad-fi
       assert(useS);
       if system.activeComponents.disgas
           % CO2 can be dissolved
           tsolve = @transport_adiOGD_simple;
       else
           tsolve = @transport_adiOG_simple;
       end
       system.nonlinear.tol = opt.transportTolerance;
       system.nonlinear.use_ecltol = false;
       
       solver = @(sol, w, bc, dT) tsolve(sol, Gt, system, bc, w, dT, fluid);
   else
       switch(lower(opt.discretization))
           case 'tpfa'
               T = computeTransmissibilityVE(Gt, rock2D);
               assert(useS);
               psolve = @(sol, w, bc) incompTPFA(sol, Gt, T, fluid,...
                                                  'wells', w, 'bc', bc);
           case 'mimetic'
               require mimetic
               S = computeMimeticIPVE(Gt, rock2D, 'Innerproduct','ip_tpf');
               if useS
                   solveMimetic = @solveIncompFlowVE;
               else
                   solveMimetic = @solveIncompFlow;
               end
               
               psolve = @(sol, w, bc) solveMimetic(sol, Gt, S, rock, fluid, ...
                                                    'bc', bc, 'wells', w);
           otherwise
               error('Unknown discretization!');
       end
       
       switch(lower(opt.transport))
           case 'implicit'
               assert(useS)
               tsolve = @(sol, w, bc, dT) implicitTransport(sol, Gt, dT, rock2D, fluid, ...
                            'wells', w, 'bc', bc, 'Trans', T, 'nltol', opt.transportTolerance);
           case 'explicit'
               if useS
                   tsolve = @(sol, w, bc, dT) explicitTransport(sol, Gt, dT,...
                                                rock, fluid,...
                                                'Trans', T, 'wells', w, 'bc', bc);
               else
                   tsolve = @(sol, w, bc, dT) explicitTransportVE(sol, Gt, dT,...
                                                rock, fluid,...
                                                'wells', w, 'bc', bc);
               end
           case 'explicit-c'
               assert(~useS)
               tsolve = @(sol, dT) mtransportWrapper(sol, Gt, dT, rock, fluid);
           otherwise
               error('Unknown transport solver');
       end
       
      solver = @(sol, w, bc, dT) transportSolveIncomp(sol, w, bc, dT, ...
                                        psolve, tsolve, Gt, fluid, useS);
   end
end


function [sol, its, conv] = transportSolveIncomp(sol, w, bc, dT, psolve, tsolve, Gt, fluid, isS)
    sol = psolve(sol, w, bc);
    if isS
        [sol, report] = tsolve(sol, w, bc, dT);
            conv.converged = report.success;
        if conv.converged
            its = report.iterations;
        else
            its = report.vasted_iterations;
        end
    else
        % We have a solver that does not report convergence.
        sol = tsolve(sol, w, bc, dT);
        conv.converged = true;
        its = 1;
    end
    
    if isS
        [s, h, hm] = normalizeValuesVE(Gt, sol, fluid); %#ok
        sol.h = h;
        sol.h_max = hm;
    else
        sol.s = height2Sat(sol, Gt, fluid);
        sol.extSat= [min(sat,sol.extSat(:,1)),max(sat,sol.extSat(:,2))];
    end

end

function T = computeTransmissibilityVE(Gt, rock2D)
       T = computeTrans(Gt, rock2D);
       T = T.*Gt.cells.H(gridCellNo(Gt));
end

function mtransportWrapper(sol, Gt, dT, rock, fluid, varargin)
      [sol.h, sol.h_max] = mtransportVE(sol, Gt, dT, rock, fluid, ...
                                          'gravity', norm(gravity));
end

%%%%%% TRANSPORT SOLVERS, FULLY IMPLICIT

function [sol, its, convergence] = transport_adiOGD_simple(sol, Gt, systemOG, bcVE_s, WVE_s, dT,fluidADI)     
    
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
    [state, its, convergence] = solvefiADI(state, dT, W, Gt, systemOG,'bc',bcVE_s);%#ok
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

function [sol, its, convergence] = transport_adiOG_simple(sol, Gt, systemOG, bcVE_s, WVE_s, dT,fluidADI)
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
    [state, its, convergence] = solvefiADI(state, dT, W, Gt, systemOG,'bc',bcVE_s);%#ok
    %[state, its] = solvefiADI(state, dT, W, Gt, systemOG,'bc',[]);%#ok
    sat=state.s(:,2);
    sol.s=sat;    
    sol.extSat= [min(sat,state.smin(:,2)),max(sat,state.smax(:,2))];
    %sol.sGmax=sol.extSat(:,2);
    sol.pressure=state.pressure;
    %
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