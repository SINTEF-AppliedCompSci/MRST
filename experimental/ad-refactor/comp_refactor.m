function comp_refactor()

    moduleCheck('ad-fi', 'ad-refactor');
    gravity on;
    
    %% Define grid, rock and temperature regime
    [Gt, rock] = makeTopSurfaceGrid([100, 1, 1],      ...  % # cells
                                    [4000, 3000 150], ...  % phys. dim.
                                    750, 0.1,         ...  % depth, porosity
                                    400 * milli*darcy);    % permeability
    ref_temp  = 273.15 + 6; % degrees kelvin
    ref_depth = 0;          % surface used as temperature reference depth
    temp_grad = 40;         % degrees per kilometer
    tinfo     = {ref_temp, ref_depth, temp_grad}; 
    rhoW      = 1020 * kilogram / meter^3; % density of brine
    
    
    %% define injection well and schedule
    tnum     = 6; %60; % total number of timesteps
    inum     = 2;%20; % number of injection steps
    tot_time = 60 * year;
    schedule = struct('W', addWell([], Gt, rock, ceil(Gt.cells.num/2), ...
                                   'type'   , 'rate'                      , ...
                                   'radius' , 0.3                         , ...
                                   'comp_i' , [0 0 1]                     , ...
                                   'val'    , 1e7 * kilo * kilogram /year , ...
                                   'name'   , 'I'), ...
                      'step', struct('val'    , diff(linspace(0, tot_time, tnum+1)), ...
                                     'control', [ones(inum,1); zeros(tnum-inum, 1)]));
    
    %% Initialize state
    % @@ NB: computation of hydrostatic pressure below must be changed if
    % reservoir is tilted
    p0 = hydrostaticPressure(Gt, rhoW, 1*atm);
    state = struct('pressure', p0, ...
                   'h', expand_var(0, Gt.cells.num), ...
                   'wellSol', struct('bhp', p0(1), 'qGs', schedule.W(1).val));
    
    %% define boundary conditions
    bc = pside([], Gt, 'LEFT',  p0(1));
    bc = pside(bc, Gt, 'RIGHT', p0(end));
    
    %% Define model
    model = fullCompressibleCO2BrineModel(Gt, rock, tinfo, 'rhoBrine', rhoW);

    %% Run schedule
    [wellSols, states] = ...
        runScheduleRefactor(state, model, schedule, 'bc', bc);

    
    %% compute caprock values and plot result
    for s = states'
        s_tmp = model.includeComputedCaprockValues(s{:});
        plot(s_tmp.h);
        pause(1);
    end
    
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

% ----------------------------------------------------------------------------
function p = hydrostaticPressure(Gt, rho_water, surface_pressure)
% ----------------------------------------------------------------------------
    p = rho_water * norm(gravity()) * Gt.cells.z + surface_pressure;
end
