function [masses, masses_0] = massTrappingDistributionVEADI(Gt, p, sG, sW, h, h_max, rock, fluidADI, trapstruct, dh, varargin)
% Compute the trapping status distribution of CO2 in each cell of a top-surface grid
%
% SYNOPSIS:
%   masses = massTrappingDistributionVEADI(Gt, p, sW, sG, h, h_max, ...
%                              rock, fluidADI, sr, sw, trapstruct)
%   masses = massTrappingDistributionVEADI(Gt, p, sW, sG, h, h_max, ...
%                              rock, fluidADI, sr, sw, trapstruct, 'rs',rs)
%
% DESCRIPTION:
%
% PARAMETERS:
%   Gt         - Top surface grid
%   p          - pressure, one value per cell of grid
%   sW         - water saturation, one value per cell of grid
%   sG         - gas saturation, one value per cell of grid
%   h          - gas height below top surface, one value per cell of grid
%   h_max      - maximum historical gas height, one value per cell of grid
%   rock       - rock parameters corresponding to 'Gt'
%   fluidADI   - ADI fluid object (used to get densities and compressibilities)
%   sr         - gas residual saturation (scalar)
%   sw         - liquid residual saturation (scalar)
%   trapstruct - trapping structure
%   dh         - subtrapping capacity (empty, or one value per grid cell of Gt)
%   varargin   - optional parameters/value pairs.  This currently only
%                includes the option 'rs', which specifies the amount of
%                dissolved CO2 (in its absence, dissolution is ignored).
%
% RETURNS:
%   masses - vector with 7 components, representing:
%            masses[1] : mass of dissolved gas, per cell
%            masses[2] : mass of gas that is both structurally and residually trapped
%            masses[3] : mass of gas that is residually (but not structurally) trapped
%            masses[4] : mass of non-trapped gas that will be residually trapped
%            masses[5] : mass of structurally trapped gas, not counting the gas that 
%                        will eventually be residually trapped
%            masses[6] : mass of subscale trapped gas (if 'dh' is nonempty)
%            masses[7] : mass of 'free' gas (i.e. not trapped in any way)
%   masses_0 (optional) - masses in terms of one value per grid cell of Gt

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

    opt.rs = 0;
    opt    = merge_options(opt, varargin{:});
    
    % Extracting relevant information from 'sol'
    sw=fluidADI.res_water;%liquid residual saturation (scalar)
    sr=fluidADI.res_gas;%gas residual saturation (scalar)
    SF = sW .* Gt.cells.H;
    SG = sG .* Gt.cells.H;
    rs = opt.rs;
    
    pvMult = 1; 
    if isfield(fluidADI, 'pvMultR')
        pvMult =  fluidADI.pvMultR(p);
    end
    pv       = rock.poro .* Gt.cells.volumes .* pvMult;%effective area
    if isfield(rock,'ntg')
       pv = rock.poro .* Gt.cells.volumes .* rock.ntg .* pvMult;
       %effective area accounting for possible net-to-gross data
    end
    rhoCO2   = fluidADI.rhoGS .* fluidADI.bG(p);
    gasPhase = sum(pv .* (rhoCO2 .* SG));
    
    if isfield(trapstruct, 'z_spill_loc')
        zt = max(trapstruct.z_spill_loc-Gt.cells.z, 0);
        zt = min(zt, Gt.cells.H);
    else
        tr   = trapstruct.trap_regions;
        tr_z = [0; trapstruct.trap_z];
        tr   = 1 + tr;
        tr(tr>numel(tr_z)) = 1;
        zt   = max(tr_z(tr) - Gt.cells.z, 0);
        zt   = min(zt, Gt.cells.H);
    end
    
    % Determine amount of CO2 in subtraps, if any
    h_sub  = zeros(Gt.cells.num, 1); % subtrapped part of 'h'
    hm_sub = zeros(Gt.cells.num, 1); % subtrapped part of 'h_max'
    if ~isempty(dh)
        h_sub  = min(dh, h);
        hm_sub = min(dh, h_max);
    end
    h_eff  = h - h_sub;
    hm_eff = h_max - hm_sub;
    
    % this requires that the fluid has a sharp interface relperm of normal type    
    hdift     = max(min(zt, hm_eff) - min(zt, h_eff),0);    % trapped part of h_max-h
    strucVol  = sum(min(zt, h_eff) .* pv .* rhoCO2);             % trapped, flowing
    plumeVol  = sum(rhoCO2 .* h_eff.* pv) - strucVol;            % non-trapped, flowing
    resStruc  = (strucVol + sum(hdift .* rhoCO2 .* pv)) * sr;      % trapped, res
    freeStruc = strucVol * (1 - sr - sw);                          % trapped, non-res
    freeRes   = plumeVol * sr;                                     % non-trapped, flowing, res
    freeMov   = plumeVol * (1 - sw - sr);                          % non-trapped, flowing, non-res
    resTrap   = sum(max(hm_eff - max(zt, h_eff),0) .* ...
                    rhoCO2 .* pv ) .* sr;                          % non-trapped, non-flowing, res
    resDis    = fluidADI.rhoGS .* sum(pv.* (rs .* fluidADI.bW(p) .* SF)); % dissolved
    subtrap   = sum((hm_sub * sr + h_sub * (1 - sr - sw)) .* pv .* rhoCO2);

    masses    = max([value(resDis), value(resStruc), value(resTrap), ...
                     value(freeRes), value(freeStruc), value(subtrap), ...
                     value(freeMov)], 0); % may be ADI variables
         
    if(abs(sum(masses(2:end))-gasPhase) > 1e-3 * gasPhase)
        disp('There is a mismatch between mass calculations');
    end
    
    if nargout > 1
        % values one per cell of grid Gt
        strucVol_0  = min(zt, h_eff) .* pv .* rhoCO2;
        plumeVol_0  = rhoCO2 .*h_eff .* pv - strucVol_0;
        resStruc_0  = (strucVol_0 + hdift.*rhoCO2.*pv) .* sr;
        freeStruc_0 = strucVol_0 .* (1-sr-sw);
        freeRes_0   = plumeVol_0 .* sr;
        freeMov_0   = plumeVol_0 .* (1-sr-sw);
        resTrap_0   = max(hm_eff - max(zt, h_eff),0) .* rhoCO2 .* pv .* sr;
        resDis_0    = fluidADI.rhoGS .* pv .* (rs .* fluidADI.bW(p) .* SF);
        subtrap_0   = (hm_sub * sr + h_sub * (1 - sr - sw)) .* pv .* rhoCO2;

        masses_0 = [{resDis_0}, {resStruc_0}, {resTrap_0}, {freeRes_0}, ...
                    {freeStruc_0}, {subtrap_0}, {freeMov_0}]; % may be ADI vars
    end
end
