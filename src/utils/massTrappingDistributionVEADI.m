function masses = massTrappingDistributionVEADI(Gt, state, rock, fluidADI, sr, sw, trapstruct, dh)
% Compute the trapping status distribution of CO2 in each cell of a top-surface grid
%
% SYNOPSIS:
%   function masses = massTrappingDistributionVEADI(Gt, sol, rock, fluidADI, sr, sw, trapstruct)
%
% DESCRIPTION:
%
% PARAMETERS:
%   Gt         - Top surface grid
%   state      - Falid 'state' object, with the additional fields 'h' and 'h_max'.
%   rock       - rock parameters corresponding to 'Gt'
%   fluidADI   - ADI fluid object (used to get densities and compressibilities)
%   sr         - gas residual saturation (scalar)
%   sw         - liquid residual saturation (scalar)
%   trapstruct - trapping structure
%   dh         - subtrapping capacity (empty, or one value per grid cell of Gt.
%
% RETURNS:
%   masses - vector with 6 components, representing:
%            masses[1] : mass of dissolved gas, per cell
%            masses[2] : mass of gas that is both structurally and residually trapped
%            masses[3] : mass of gas that is residually (but not structurally) trapped
%            masses[4] : mass of non-trapped gas that will be residually trapped
%            masses[5] : mass of structurally trapped gas, not counting the gas that 
%                        will eventually be residually trapped
%            masses[6] : mass of 'free' gas (i.e. not trapped in any way)
%            masses[7] : mass of subscale trapped gas (if 'dh' is nonempty)

    % Extracting relevant information from 'sol'
    p = state.pressure;
    sF = state.s(:,1);
    sG = state.s(:,2);
    SF = sF .* Gt.cells.H;
    SG = sG .* Gt.cells.H;
    rs = 0;
    if isfield(state, 'rs')
        rs = state.rs;
    end
    
    pvMult = 1; 
    if isfield(fluidADI, 'pvMultR')
        pvMult =  fluidADI.pvMultR(p);
    end
    pv       = rock.poro .* Gt.cells.volumes .* pvMult;
    rhoCO2   = fluidADI.rhoG .* fluidADI.bG(p);
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
        h_sub  = min(dh, state.h);
        hm_sub = min(dh, state.h_max);
    end
    h_eff  = state.h - h_sub;
    hm_eff = state.h_max - hm_sub;
    
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
    resDis    = fluidADI.rhoG .* sum(pv .* (rs .* fluidADI.bO(p) .* SF)); % dissolved
    subtrap   = sum((hm_sub * sr + h_sub * (1 - sr - sw)) .* pv .* rhoCO2);

    masses    = max([resDis, resStruc, resTrap, freeRes, freeStruc,  freeMov, subtrap], 0);

    if(abs(sum(masses(2:end))-gasPhase) > 1e-3 * gasPhase)
        disp('There is a mismatch between mass calculations');
    end
end

    
