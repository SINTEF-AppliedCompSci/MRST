function very_simple_example(varargin)

   gravity on;
   
   Gt = setup_grid();
   rock = setup_rock();
   fluid = setup_fluid();
   state = setup_state();
   model = CO2VEBlackOilTypeModel(Gt, rock, fluid);
   
   schedule = setup_schedule();
   [wellSols, states] = simulateScheduleAD(state, model, schedule);
   
end
% ============================================================================



