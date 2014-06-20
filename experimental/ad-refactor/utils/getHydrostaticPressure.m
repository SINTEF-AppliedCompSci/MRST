function p = getHydrostaticPressure(Gt, rho_water, surface_pressure, slope, ...
                                 slopedir, h)
% ----------------------------------------------------------------------------
    
    depth = computeRealDepth(Gt, slope, slopedir, h);
    
    p = rho_water * norm(gravity()) * depth + surface_pressure;
end
