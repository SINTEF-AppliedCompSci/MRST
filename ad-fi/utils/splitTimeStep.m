function schedule = splitTimeStep(schedule, step)
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
