function comp_refactor()

    moduleCheck('ad-fi');
    gravity on;
    
    %% Define grid, rock and temperature regime
    [Gt, rock] = makeTopSurfaceGrid([100, 1, 1],      ...  % # cells
                                    [4000, 3000 150], ...  % phys. dim.
                                    750, 0.1,         ...  % depth, porosity
                                    400 * milli*darcy);    % permeability
    ref_temp  = 273.15 + 6; % degrees kelvin
    ref_depth = 0;          % surface used as temperature reference depth
    temp_grad = 40;         % degrees per kilometer
    tinfo     = {ref_temp, ref_dept, temp_grad}; 
    rhoW      = 1020 kilogram / meter^3; % density of brine
        
    %% define injection well and schedule
    schedule = ;
    
    %% Initialize state
    state = struct('pressure', hydrostaticPressure(Gt, rhoW, 1*atm), ...
                   'h'       , expand_var(0, Gt.cells.num));
    
    %% Define model
    model = fullCompressibleCO2BrineModel(Gt, tinfo, 'rhoBrine', rhoW);
    
    %% Define solver
    nlsolve = adaptiveNonlinearSolver('linearSolver', CPRSolverAD());
    
    %% Run schedule
    [wellSols, states] = ...
        runScheduleRefactor(state, model, schedule, 'nonlinearSolver', nlsolve);
end

% ----------------------------------------------------------------------------
function [Gt, rock] = makeTopSurfaceGrid(dim, lengths, depth, poro, perm)
% ----------------------------------------------------------------------------
    Gt   = topSurfaceGrid(computeGeometry(setTopDepth(cartGrid(dim, lengths), depth)));
    rock = averageRock(struct('perm', expand_var(perm, prod(dim)), ...
                              'poro', expand_var(poro, prod(dim))), Gt);
end

% ----------------------------------------------------------------------------
function G = setTopDepth(G, depth)
% ----------------------------------------------------------------------------
    mindepth = min(G.nodes.coords(:,3));
    G.nodes.coords(:,3) = G.nodes.coords(:,3) + (depth - mindepth);
end
        