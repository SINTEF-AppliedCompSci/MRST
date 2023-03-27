function [Gt, wellSols, states] = very_simple_example(varargin)

   gravity on;
   
   [G, Gt]            = setup_grid             ();
   rock               = setup_rock             (Gt);
   fluid              = setup_fluid            (Gt, rock);
   [state, bc]        = setup_state_and_bc     (Gt);
   model              = CO2VEBlackOilTypeModel (Gt, rock, fluid);
   schedule           = setup_schedule         (G, Gt, rock, bc);
   [wellSols, states] = simulateScheduleAD     (state, model, schedule);
   
end
% ============================================================================

function [G, Gt] = setup_grid()

   G = cartGrid([5 5 1], [10 10 0.1] * kilo * meter);
   G = computeGeometry(G);
   Gt = topSurfaceGrid(G);
end

function rock = setup_rock(Gt)
   
   rock = struct('perm', 1000 * milli * darcy * ones(Gt.cells.num, 1), ...
                 'poro', 0.2 * ones(Gt.cells.num, 1));
   rock = averageRock(rock, Gt);
end

function fluid = setup_fluid(Gt, rock)
   
   fluid = makeVEFluid(Gt, rock, 'simple'); 
   
end

function [state, bc] = setup_state_and_bc(Gt)
   
   pressure = 100 * barsa;

   bc = pside([], Gt, 'XMin', pressure);
   bc = pside(bc, Gt, 'XMax', pressure);
   bc = pside(bc, Gt, 'YMin', pressure);
   bc = pside(bc, Gt, 'Ymax', pressure);

   bc.sat = [ones(size(bc.value)), zeros(size(bc.value))]; % @@ should not need to
                                                           % be set manually
   state.pressure = pressure * ones(Gt.cells.num, 1);
   state.s = [ones(Gt.cells.num, 1), zeros(Gt.cells.num, 1)];
   state.sGmax = state.s(:,2);
   
end

function schedule = setup_schedule(G, Gt, rock, bc)

   cell_ix = 13;
   W = addWell([], G, rock, cell_ix, ...
               'Type', 'bhp', ...
               'val', 120 * barsa, ...%120 * barsa, ...
               'Radius', 0.3, ...
               'comp_i', [0 1], ...
               'name', 'I');

   W.sign = 1; % @@ Should not have to be manually handled

   W = convertwellsVE(W, G, Gt, rock);
   
   schedule = struct('control', struct('W', W, 'bc', bc), ...
                     'step', struct('control', ones(10, 1), ...
                                    'val', year * ones(10,1)));
end
