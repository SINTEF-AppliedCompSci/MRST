function new_schedule = refineSchedule(init_time_for_new_time_steps, ...
                                   new_time_steps, ...
                                   schedule)

   % Compute a finer schedule, including new time steps but preserving the time steps of the original
   % schedule
   %
   % SYNOPSIS: 
   %  new_schedule = refineSchedule(init_time_for_new_time_steps, new_time_steps, schedule)
   %
   %
   % PARAMETERS:
   %
   %   init_time_for_new_time_steps  - Time where the sequence of new time steps will be added.
   %   
   %   new_time_steps                - Sequence of time steps to be added.
   %
   %   schedule                      - Input schedule which will be refined
   %

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
   
   ref_steps = schedule.step.val;
   ref_control = schedule.step.control;

   n1 = numel(new_time_steps) + 1;
   n2 = numel(ref_steps) + 1;
   cum_time = [[init_time_for_new_time_steps + [0; cumsum(new_time_steps)], ones(n1, 1)]; [init_time_for_new_time_steps + [0; ...
                       cumsum(ref_steps)], 2*ones(n2, 1)]];
   cum_time  = sortrows(cum_time, 1);
   n = size(cum_time, 1);


   % remove elements which are identical (close) to reference schedule
   remove_ind = [];
   epsi     = 1; % 1 second
   for i = 1 : n 
      if cum_time(i, 2) == 2
         t = cum_time(i, 1);
         if i > 1
            tm = cum_time(i - 1, :);
            if tm(2) == 1 && abs(tm(1) - t) < epsi
               remove_ind = [remove_ind; i - 1];
            end
         end
         if i < n
            tp = cum_time(i + 1, :);
            if tp(2) == 1 && abs(tp(1) - t) < epsi
               remove_ind = [remove_ind; i + 1];
            end
         end
      end
   end
   ind = true(n, 1);
   ind(remove_ind) = false;
   cum_time = cum_time(ind, 1);
            
   steps = diff(cum_time);
   
   t        = init_time_for_new_time_steps;
   ref_time = init_time_for_new_time_steps + ref_steps(1);
   ref_ind  = 1;
   N        = numel(steps);
   control  = zeros(N, 1);

   for ind = 1 : N
      t = t + steps(ind);
      while (t - ref_time > epsi)
         ref_ind = ref_ind + 1;
         ref_time = ref_time + ref_steps(ref_ind);
      end
      control(ind) = ref_control(ref_ind);
   end
   
   if isfield(schedule, 'control')
      new_schedule.control = schedule.control;
   end
   if isfield(schedule, 'W')
      new_schedule.W = schedule.W;
   end
   if isfield(schedule, 'type')
      new_schedule.type = schedule.type;
   end
   new_schedule.step.val = steps;
   new_schedule.step.control = control;
   
   new_schedule.step.repStep = true(numel(new_schedule.step.val), 1);
   
end
