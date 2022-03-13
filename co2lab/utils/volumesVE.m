function vol = volumesVE(G, sol, rock, fluid, ts)
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
%        moved out, defined as the difference between h_max and h
%      - the free residual volume defined as the height of the CO2
%        colum in each cell multiplied by the pore volume of the cell
%        and the residual CO2 saturation
%      - free volume defined as the height of CO2 column in each cell
%        multiplied by the pore volume of the cell and the CO2
%        saturation minus the residual CO2 saturation
%   If a trapping structure is provided, the vector consists of five
%   entries:
%      - residual volumes of CO2 confined to structural traps
%      - residual volume of CO2 left in cells the CO2 plume has
%        moved out of
%      - fraction of the free volume that will remain as residual CO2
%        when the plume moves away from its current location
%      - movable volumes of CO2 confined in structural traps
%      - fraction of the CO2 plume that is free to leave the current
%        position (i.e., will not remain as residually trapped)
%  
%        The five entries of the vector can be summarized using the below table:
%
%                                        -------- entry of 'vol' -------
%        | Zone/type                     | 1   | 2   | 3   | 4   | 5   |
%        |-------------------------------+-----+-----+-----+-----+-----|
%        | In trap                       | yes | no  | no  | yes | no  |
%        | Outside trap                  | no  | yes | yes | no  | yes |
%        | Residual                      | yes | yes | yes | no  | no  |
%        | In presence of free flow (<h) | y/n | no  | yes | yes | yes |
%        |-------------------------------+-----+-----+-----+-----+-----|
%

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

    pv = rock.poro.*G.cells.volumes;
    if nargin==4
       trapped  = sum((sol.h_max-sol.h).*pv)*fluid.res_gas;
       plumeVol = sum(sol.h.*pv);
       freeRes  = plumeVol*fluid.res_gas;
       free     = plumeVol*(1-fluid.res_water-fluid.res_gas);
       vol      = [trapped freeRes free];
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
       resStruc  = (strucVol + sum(hdift.*pv)) * fluid.res_gas;
       freeStruc = strucVol*(1 - fluid.res_gas - fluid.res_water);
       freeRes   = plumeVol*fluid.res_gas;
       freeMov   = plumeVol*(1-fluid.res_water - fluid.res_gas);
       resTrap   = sum(max(sol.h_max - max(zt, sol.h),0).*pv ).*fluid.res_gas;
       vol       = max([resStruc, resTrap, freeRes, freeStruc, freeMov],0);
    end
end
