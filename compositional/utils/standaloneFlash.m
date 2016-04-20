function [L, x, y, Z_L, Z_V] = standaloneFlash(p, T, z, EOSModel)
% Utility for flashing without explicitly forming a state
    solver = NonLinearSolver();
    solver.maxTimestepCuts = 0;
    solver.maxIterations = 200;
    
    state = struct();

    state.pressure = p;
    state.T = T;
    state.components = z;
    
    state = solver.solveTimestep(state, 1, EOSModel);
    
    L = state.L;
    x = state.x;
    y = state.y;
    Z_L = state.Z_L;
    Z_V = state.Z_V;
end