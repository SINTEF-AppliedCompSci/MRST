function reports = makeReports(Gt, states, rock, fluid, schedule, ...
                               residual, traps, dh)
% residual on form: [sw, sr]
   
   tot_inj = 0;
   for i = 1:numel(states)

      [h, h_max] = computePlumeHeight(Gt, states{i}, residual(1), residual(2));

      if i == 1
         reports(i).t         = 0;
         reports(i).W         = [];
      else
         reports(i).t = sum(schedule.step.val(1:i-1));
         reports(i).W = schedule.control(schedule.step.control(i-1)).W;
         
         assert(all(cellfun(@(x) strcmpi(x, 'rate'), {reports(i).W.type})));
         tot_inj = tot_inj + sum([reports(i).W.val]) * schedule.step.val(i-1) * fluid.rhoGS;
      end
      
      reports(i).sol       = states{i};
      reports(i).sol.h     = h;
      reports(i).sol.h_max = h_max;
      
      reports(i).masses    = massTrappingDistributionVEADI(Gt          , ...
                                                        reports(i).sol , ...
                                                        rock           , ...
                                                        fluid          , ...
                                                        traps          , ...
                                                        dh);
      leaked = tot_inj - sum(reports(i).masses);
      reports(i).masses = [reports(i).masses, leaked];
   end
   
end

