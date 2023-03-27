function masses = massesVE(G, sol, rock, fluid, ts)
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
   
    pv = rock.poro.*G.cells.volumes;
    if nargin==4
       trapped = sum((sol.h_max-sol.h).*pv)*fluid.sr;
       free    = sum(sol.h.*pv)*(1-fluid.sw);
       vol     = [trapped free];
    elseif nargin==5
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
       hdift     = max(min(zt, sol.h_max) - min(zt, sol.h),0);
       strucVol  = sum(min(zt, sol.h).*pv);
       plumeVol  = sum(sol.h.* pv) - strucVol;
       resStruc  = (strucVol + sum(hdift.*pv)) * fluid.sr;
       freeStruc = strucVol*(1 - fluid.sr - fluid.sw);
       freeRes   = plumeVol*fluid.sr;
       freeMov   = plumeVol*(1-fluid.sw - fluid.sr);
       resTrap   = sum(max(sol.h_max - max(zt, sol.h),0).*pv ).*fluid.sr;
       resDis    = 0;
       masses       = fluid.rho(1)*max([resDis, resStruc, resTrap, freeRes, freeStruc, freeMov,],0);
    end
end
