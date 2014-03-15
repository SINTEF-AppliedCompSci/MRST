function masses = massTrappingDistributionVEADI(Gt, sol, rock, fluidADI, sr, sw, trapstruct)
% Compute the trapping status distribution of CO2 in each cell of a top-surface grid
%
% SYNOPSIS:
%   function masses = massTrappingDistributionVEADI(Gt, sol, rock, fluidADI, sr, sw, trapstruct)
%
% DESCRIPTION:
%
% PARAMETERS:
%   Gt         - Top surface grid
%   sol        - Solution structure, containing a valid 'state' object, as
%                well as 'h' and 'h_max' fields.
%   rock       - rock parameters corresponding to 'Gt'
%   fluidADI   - ADI fluid object (used to get densities and compressibilities)
%   sr         - gas residual saturation (scalar)
%   sw         - liquid residual saturation (scalar)
%   trapstruct - trapping structure
%
% RETURNS:
%   masses - vector with 6 components, representing:
%            masses[1] : mass of dissolved gas, per cell
%            masses[2] : mass of gas that is both structurally and residually trapped
%            masses[3] : mass of gas that is residually (but not structurally) trapped
%            masses[4] : mass of non-trapped gas that will be residually trapped
%            masses[5] : mass of structurally trapped gas, not counting the gas that 
%                        will eventually be residually trapped
%            masses[6] : mass of structurally (but not residually) trapped gas
%            masses[7] : mass of 'free' gas (i.e. not trapped in any way)

    % Extracting relevant information from 'sol'
    p = sol.state.pressure;
    sF = sol.state.s(:,1);
    sG = sol.state.s(:,2);
    SF = sF .* Gt.cells.H;
    SG = sG .* Gt.cells.H;
    rs = 0;
    if isfield(sol.state, 'rs')
        rs = sol.state.rs;
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
    
    % this depend that the fluid has a sharp interface relperm of normal
    % type
    hdift     = max(min(zt, sol.h_max) - min(zt, sol.h),0);
    strucVol  = sum(min(zt, sol.h) .* pv .* rhoCO2);
    plumeVol  = sum(rhoCO2 .* sol.h.* pv) - strucVol;
    resStruc  = (strucVol + sum(hdift .* rhoCO2 .* pv)) * sr;
    freeStruc = strucVol * (1 - sr - sw);
    freeRes   = plumeVol * sr;
    freeMov   = plumeVol * (1 - sw - sr);
    resTrap   = sum(max(sol.h_max - max(zt, sol.h),0) .* rhoCO2 .* pv ) .* sr;
    resDis    = fluidADI.rhoG .* sum(pv .* (rs .* fluidADI.bO(p) .* SF)); 

    masses    = max([resDis, resStruc, resTrap, freeRes, freeStruc,  freeMov],0);

    if(abs(sum(masses(2:end))-gasPhase) > 1e-3 * gasPhase)
        disp('There is a mismatch between mass calculations')
    end
end

    
