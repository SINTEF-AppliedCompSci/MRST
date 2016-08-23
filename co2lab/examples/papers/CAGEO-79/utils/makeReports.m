function reports = makeReports(Gt, states, rock, fluid, schedule, ...
                               residual, traps, dh)
   % This function does intermediate processing of simulation data for
   % further visualization by 'selectedResultsMultiplot'.

   % residual on form: [s_water, s_co2]
   % states on form:   {initState, states{:}} to include initial state

   assert( numel(states) == numel(schedule.step.val)+1 , ...
       'Ensure the initial state has been included in the varargin ''states''.')
   
   tot_inj = 0;
   for i = 1:numel(states)

      [h, h_max] = computePlumeHeight(Gt, states{i}, residual(1), residual(2));

      if i == 1
         reports(i).t         = 0; %#ok
         reports(i).W         = []; %#ok
      else
         reports(i).t = sum(schedule.step.val(1:i-1)); %#ok
         reports(i).W = schedule.control(schedule.step.control(i-1)).W; %#ok
         
         assert(all(cellfun(@(x) strcmpi(x, 'rate'), {reports(i).W.type})));
         tot_inj = tot_inj + sum([reports(i).W.val]) * schedule.step.val(i-1) * fluid.rhoGS;
      end
      
      reports(i).sol       = states{i}; %#ok
      reports(i).sol.h     = h; %#ok
      reports(i).sol.h_max = h_max; %#ok
      
      reports(i).masses    = massTrappingDistributionVEADI(Gt                   , ...
                                                        reports(i).sol.pressure , ...
                                                        reports(i).sol.s(:,2)   , ...
                                                        reports(i).sol.s(:,1)   , ...
                                                        reports(i).sol.h        , ...
                                                        reports(i).sol.h_max    , ...
                                                        rock                    , ...
                                                        fluid                   , ...
                                                        traps                   , ...
                                                        dh);%#ok
      
      leaked = tot_inj - sum(reports(i).masses);
      reports(i).masses = [reports(i).masses, leaked]; %#ok
   end
   
end

