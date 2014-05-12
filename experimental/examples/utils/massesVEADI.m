function masses = massesVEADI(G, sol, rock, fluidADI, fluid,ts)
    error(['massesVEADI is deprecated.  Use phaseMassesVEADI or ' ...
           'MassTrappingDistributionVEADI instead.\n']);
    % SYNOPSIS:
    %   vol = volumesVE(G, sol, rock, fluid)
    %   vol = volumesVE(G, sol, rock, fluid, ts)
    %   
    % PARAMETERS:
    %   G      - 2D top surface grid used for VE-simulations
    %   sol    - Solution state as defined by initResSolVE
    %   rock   - rock for the top surface grid
    %   fluid  - fluid object
    %   ts     - trapping structure object, typically returned by a call to
    %            'findTrappingStructure'
    % 
    % RETURNS:
    %   a vector with trapped and free volumes. If no trapping structure is
    %   provided, the vector consists of two entries:
    %      - the residual CO2 saturation in regions where the CO2 plume has
    %        moved out again defined as the difference between h_max and h
    %      - free volume defined as the height of CO2 column in each cell
    %        multiplied by the pore volume of the cell and the CO2
    %        saturation
    %   If a trapping structure is provided, the vector consists of four
    %   entries:
    %      - residual volumes of CO2 confined to structural traps
    %      - residual volume of CO2 left in cells the CO2 plume has
    %        moved out again
    %      - fraction of the free volume that will remained as residual CO2
    %        when the plume moves out of cells
    %      - movable volumes of CO2 confined to structural traps
    %      - fraction of the free volume that can leave the cells
   pvMult = 1; 
   if isfield(fluidADI, 'pvMultR')
    pvMult =  fluidADI.pvMultR(sol.state.pressure);
   end
   if(isfield(sol.state,'sGmax'))
    pcOG  = fluidADI.pcOG(sol.state.s(:,2),sol.state.pressure,'sGmax',sol.state.sGmax);
   else
     pcOG  = fluidADI.pcOG(sol.state.s(:,2),sol.state.pressure,'sGmax',sol.state.smax(:,2));  
   end
   rhoCO2=fluidADI.rhoG.*fluidADI.bG(sol.state.pressure);%+pcOG);
    pv = rock.poro.*G.cells.volumes.*pvMult;
    gasPhase =  sum(pv.*(rhoCO2.*G.cells.H.*(sol.state.s(:,2))));
    %%
    %gasPhase_h = sum(pv.*(rhoCO2.*(sol.h).*(1-fluid.sw)));
    %gasPhase_h_max = sum(pv.*(rhoCO2.*(sol.h_max -sol.h).*(fluid.sr)));
    
    if nargin==5
        if(false)
        trapped = sum((sol.h_max-sol.h).*pv)*fluid.sr;
        free    = sum(sol.h.*pv)*(1-fluid.sw);
        vol     = [trapped free];
        else
          rhoW=fluidADI.rhoO.*fluidADI.bO(sol.state.pressure);%,sol.state.rs,sol.state.s(:,2)>0);
          %gasPhase =  sum(pv.*(fluidADI.bG(sol.state.pressure).*(sol.state.s(:,2)))); 
          watPhase =  sum(pv.*(rhoW.*G.cells.H.*(sol.state.s(:,1))));  
          if(isfield(sol.state,'rs'))
            resDis    = fluidADI.rhoG.*sum(pv.*(sol.state.rs.*fluidADI.bO(sol.state.pressure)...
                .*G.cells.H.*(sol.state.s(:,1))));   
          else
             resDis=0;
          end
          masses = [gasPhase,watPhase,resDis];  
        end
        
    elseif nargin==6
       if isfield(ts, 'z_spill_loc')
           zt = max(ts.z_spill_loc-G.cells.z, 0);
           zt = min(zt, G.cells.H);
       else
           tr   = ts.trap_regions;
           tr_z = [0; ts.trap_z];
           tr   = 1 + tr;
           tr(tr>numel(tr_z)) = 1;
           zt   = max(tr_z(tr) - G.cells.z, 0);
           zt   = min(zt, G.cells.H);
       end
       
       
       % this depend that the fluid has a sharp interface relperm of normal
       % type
       hdift     = max(min(zt, sol.h_max) - min(zt, sol.h),0);
       strucVol  = sum(min(zt, sol.h).*pv.*rhoCO2);
       plumeVol  = sum(rhoCO2.*sol.h.* pv) - strucVol;
       resStruc  = (strucVol + sum(hdift.*rhoCO2.*pv)) * fluid.sr;
       freeStruc = strucVol*(1 - fluid.sr - fluid.sw);
       freeRes   = plumeVol*fluid.sr;
       freeMov   = plumeVol*(1-fluid.sw - fluid.sr);
       resTrap   = sum(max(sol.h_max - max(zt, sol.h),0).*rhoCO2.*pv ).*fluid.sr;
       
       resDis    = fluidADI.rhoG.*sum(pv.*(sol.state.rs.*fluidADI.bO(sol.state.pressure).*G.cells.H.*(sol.state.s(:,1)))); 
       
       masses       = max([resDis, resStruc, resTrap, freeRes, freeStruc, freeMov],0);
       if(abs(sum(masses(2:end))-gasPhase)>1e-3*gasPhase)
           disp('There is a miss match between mass calculations')
       end
    end
end
