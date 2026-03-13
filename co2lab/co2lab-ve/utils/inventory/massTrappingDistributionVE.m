function masses = massTrappingDistributionVE(S, Smax, pressure, rs, Gt, fluidVE, rockVE, trapstruct, dh)
% Compute the trapping distribution of CO2 in each cell of a top-surface grid
%
% SYNOPSIS:
%   function masses = massTrappingDistributionVE(S, Smax, pressure, rs, Gt, ...
%                                                fluidVE, rockVE, trapstruct, dh)
%
% DESCRIPTION:
%
% PARAMETERS:
%   S          - current CO2 saturation per cell
%   Smax       - historically maximum OC2 saturation per cell
%   pressure   - cell-wise pressure
%   rs         - fraction of dissolved CO2 in brine per cell
%   Gt         - top surface grid
%   fluidVE    - VE fluid object
%   rockVE     - VE rock object
%   trapstruct - trapping structure
%   dh         - subtrapping capacity (empty, on one value per grid cell)
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

%{
Copyright 2009-2025 SINTEF Digital, Mathematics & Cybernetics.

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

% NB: Key approximation done: for capillary fringe models, it is assumed that
% the mobile and immobile part of the saturations are uniformly distributed
% in z in [0, h] and [0,hmax], respectively.  It is also assumed that the
% enpoint scaling model is used for hysteresis, in order for the computation
% of S_immob to make sense below.

    if isempty(dh)
        dh = 0;
    end
        
    % pick relevant variables from objects
    rw = fluidVE.res_water;
    rg = fluidVE.res_gas;
    H = Gt.cells.H;
    zt = extract_spillpoint_vals(trapstruct, Gt);
    
    % compute pore volumes and densities
    pvMult = 1;
    if isfield(fluidVE, 'pvMultR')
        pvMult = fluidVE.pvMultR(pressure);
    end
    pv = rockVE.poro .* Gt.cells.volumes .* H .* pvMult; % 3D pore volumes
    if isfield(rockVE, 'ntg')
        pv = pv .* rockVE.ntg;
    end
    rhoG = fluidVE.rhoGS .* fluidVE.bG(pressure);
    mfac = pv .* rhoG;
    
    total_undissolved_mass = sum(mfac .* S); % used for sanity check
    
    % compute interface positions
    [h, h_max] = compute_plume_height(S, Smax, pressure, Gt, fluidVE, rockVE);
    h = h + eps;
    h_max = h_max + eps; % avoid division by zero in formulas below
    
    % determine free (mobile) and immobile part of current saturation
    S_immob = Smax * rg / (1-rw); % remaining saturation when mobile saturation reaches 0 (see formula in free_sg.m)
    S_mob = S - S_immob; %free_sg(S, Smax, rw, rg);
    
    % determine subscale trapping
    S_immob_sub = S_immob .* min(dh, h_max) ./ h_max;
    S_mobile_sub = S_mob .* min(dh, h) ./ h;
    
    % determine mobile and immobilized CO2 inside structural traps
    S_immob_struct = S_immob .* max(min(h_max, zt) - dh, 0) ./ h_max;
    S_mobile_struct = S_mob .* max(min(h, zt) - dh, 0) ./ h;
    
    % determine mobile and immobilized CO2 outside structural traps
    S_immob_free = S_immob - S_immob_struct - S_immob_sub;
    S_mobile_free = S_mob - S_mobile_struct - S_mobile_sub;
    
    % determine "residual trapping" (immobilized CO2 outside the flowing 
    % plume, and outside of traps).  This is a subset of 'S_immob_free'.
    S_res = S_immob .* max((h_max - max(h, zt)), 0) ./ h_max;
    
    % --- computing masses to go into reporting ---
    
    % determine dissolved CO2
    disgas_mass = sum(fluidVE.rhoGS .* fluidVE.bW(pressure) .* pv .* rs .* (1-S));
    
    % immobilized structural trapped CO2
    immob_struct_mass = sum(mfac .* S_immob_struct);
    
    % residual CO2
    residual_mass = sum(mfac .* S_res);
    
    % immobilized, free CO2 (other than residual, i.e. inside plume)
    immob_free_mass = sum(mfac .* S_immob_free) - residual_mass;
    
    % mobile structural trapped CO2
    mobile_struct_mass = sum(mfac .* S_mobile_struct);
    
    % CO2 in subtraps
    subtrap_mass = sum(mfac .* (S_immob_sub + S_mobile_sub));
    
    % mobile, free-flowing CO2
    mobile_free_mass = sum(mfac .* S_mobile_free);
    
    % assemble return variable
    masses = max([value(disgas_mass), ...         
                  value(immob_struct_mass), ...   
                  value(residual_mass), ...
                  value(immob_free_mass), ... 
                  value(mobile_struct_mass), ...  
                  value(subtrap_mass), ...        
                  value(mobile_free_mass)], 0);   
                  
    
    % sanity check
    if(abs(sum(masses(2:end))-total_undissolved_mass) > 1e-3 * total_undissolved_mass)
        abs(sum(masses(2:end))-total_undissolved_mass)
        abs(sum(masses(2:end))-total_undissolved_mass) - 1e-3 * total_undissolved_mass
        disp('There is a mismatch between mass calculations');
    end

end

% ----------------------------------------------------------------------------
function [h, h_max] = compute_plume_height(S, Smax, pressure, Gt, fluidVE, rockVE)
    
    H = Gt.cells.H;

    if is_capillary_fringe_model(fluidVE)
        % capillary fringe model
        [h, h_max] = upscaledSat2height(S, Smax, Gt, ...
                                        'pcWG', fluidVE.pcWG, ...
                                        'rhoW', fluidVE.rhoW, ...
                                        'rhoG', fluidVE.rhoG, ...
                                        'p', pressure);
        h = min(h, Gt.cells.H);
        h_max = min(h_max, Gt.cells.H);
    else
        % sharp interface model
        rg = fluidVE.res_gas;
        rw = fluidVE.res_water;
        poro = [];
        if strcmpi(fluidVE.relperm_model, 'sharp_interface_integrated')
            poro = rockVE.parent.poro;
        end
        [h, h_max] = upscaledSat2height(S, Smax, Gt, ...
                                        'resSat', [rw, rg], ...
                                        'poro', poro);
    end
end

% ----------------------------------------------------------------------------
function res = is_capillary_fringe_model(fluidVE)
    % check if this is a fluidVE model that incorporates a capillary fringe
    res = strcmpi(fluidVE.relperm_model, 'P-scaled table') || ...
          strcmpi(fluidVE.relperm_model, 'P-K-scaled table') || ...
          strcmpi(fluidVE.relperm_model, 'S table');
end

% ----------------------------------------------------------------------------
function zt = extract_spillpoint_vals(trapstruct, Gt)
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
end
