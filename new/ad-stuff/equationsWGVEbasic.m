function [problem, state] = equationsWGVEbasic(state0, state, model, dt, drivingForces, varargin)
   
   opt = struct('reverseMode', false, ...
                'resOnly', false);
   
   assert(isempty(drivingForces.src)); % unsupported
   W  = drivingForces.Wells;
   bc = drivingForces.bc; 
   s  = model.operators; 
   Gt = model.G; 
   f  = model.fluid; 
   
   
   [p, sat, wellSol] = model.getProps(state, 'pressure', 'saturation', 'wellsol');

   
   %% @@ ??? where is sGmax updated? 
end
