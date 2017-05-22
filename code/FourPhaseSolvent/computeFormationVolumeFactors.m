function [bW_eff, bO_eff, bG_eff, bS_eff] = computeFormationVolumeFactors(fluid, p, rhoO_eff, rhoG_eff, rhoS_eff)
    % Calculate effective (inverse) formation volume factors due to new densities
    
    % Immiscible formation volume factors b = rho/rhoS
    bO_i = fluid.bO(p);
    bG_i = fluid.bG(p);
    bS_i = fluid.bS(p);
    
    % Miscible formation volume factors b = rho_eff/rhoS
    bO_m = rhoO_eff./fluid.rhoOS;
    bG_m = rhoG_eff./fluid.rhoGS;
    bS_m = rhoS_eff./fluid.rhoSS;
    
    % Effective formation volume factors are interpolated using a
    % pressure-dependent miscibility funciton
    Mp = fluid.Mpres(p);
    bW_eff = fluid.bW(p);
    bO_eff = bO_m.*Mp + bO_i.*(1-Mp);
    bG_eff = bG_m.*Mp + bG_i.*(1-Mp);
    bS_eff = bS_m.*Mp + bS_i.*(1-Mp);
    
end