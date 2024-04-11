function stateVE = state2VE(state3D, Gt, fluid, res_wat, res_gas, poro3D)

    [S, Smax] = finescale2upscaledSat(state3D.s(:,2), res_wat, res_gas, Gt, poro3D);
    
    stateVE.s = [1-S, S];
    stateVE.sGmax = Smax;
    stateVE.pressure = finescale2upscaledPressure(state3D.pressure, Gt, fluid);
    
end
