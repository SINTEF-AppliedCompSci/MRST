function schedule = splitTimeStep(schedule, step)
%Undocumented Utility Function

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

   time_step = schedule.step.val(step);
   control = schedule.step.control(step);
   if step == 1
      step_val1 = [];
      step_control1 = [];
   else
      step_val1 = schedule.step.val(1:step - 1);
      step_control1 = schedule.step.control(1:step - 1);
   end
   if step == numel(schedule.step.val)
      step_val2 = [];
      step_control2 = [];
   else
      step_val2 = schedule.step.val(step + 1:end);
      step_control2 = schedule.step.control(step + 1:end);
   end
   schedule.step.val =  [step_val1; time_step/2; time_step/2; step_val2];
   schedule.step.control =  [step_control1; control; control; step_control2];
   schedule.step.repStep = true(numel(schedule.step.val), 1);
   
end
