function [bW_eff, bO_eff, bG_eff, bS_eff] = computeFormationVolumeFactors(fluid, p, rhoO_eff, rhoG_eff, rhoS_eff)
    % Calculate effective (inverse) formation volume factors due to new densities
    
    % New formation volume factors b = rho_eff/rhoS
    bW_eff = fluid.bW(p);
    bO_eff = rhoO_eff./fluid.rhoOS;
    bG_eff = rhoG_eff./fluid.rhoGS;
    bS_eff = rhoS_eff./fluid.rhoSS;
    
end