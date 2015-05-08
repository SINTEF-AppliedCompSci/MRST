function p = computeHydrostaticPressure(Gt, rho_water, surface_pressure)
    
    p = rho_water * norm(gravity()) * Gt.cells.z + surface_pressure;
end
