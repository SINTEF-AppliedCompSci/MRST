function new_schedule = refineSchedule(init_time_for_new_time_steps, ...
                                   new_time_steps, ...
                                   schedule, ...
                                   varargin)

   % Compute a finer schedule, including new time steps but preserving the time steps of the original
   % schedule
   %
   % SYNOPSIS: 
   %  new_schedule = refineSchedule(init_time_for_new_time_steps, new_time_steps, ref_init_time, schedule)
   %  new_schedule = refineSchedule(init_time_for_new_time_steps, new_time_steps, ref_init_time,
   %  schedule, 'pn', pv)
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
   %   'pn'/pv                       - List of 'key'/value pairs defining optional parameters.  The
   %                                   supported options are:
   %           
   %   ref_init_time                 - Starting time for the schedule to de modified (if
   %                                   different than zero. Useful when dealing with 
   %                                   restarted schedules.)
   
   opt = struct('ref_init_time', 0);
   opt = merge_options(opt, varargin{:});
   ref_init_time = opt.ref_init_time;
   
   ref_steps = schedule.step.val;
   ref_control = schedule.step.control;

   n1 = numel(new_time_steps) + 1;
   n2 = numel(ref_steps) + 1;
   cum_time = [[init_time_for_new_time_steps + [0; cumsum(new_time_steps)], ones(n1, 1)]; [ref_init_time + [0; ...
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
   
   t        = ref_init_time;
   ref_time = ref_init_time + ref_steps(1);
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
