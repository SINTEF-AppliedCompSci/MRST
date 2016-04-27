function obj = account_well_cost_rerun(Gt, fluid, wellSols, schedule, initial_cost, operation_rate, tax_credit_rate, alpha, nonlinear_well_cost, varargin )
%function obj = account_well_cost(model, wellSols, schedule, initial_cost, operation_rate, tax_credit_rate, alpha, varargin )
% operation_rate is in units of USD per tonne CO2
% tax_credit_rate is in units of USD per tonne CO2
% initial_cost is in units of USD
% obj{steps} returned is in units of USD

   opt.ComputePartials = false;
   opt.tStep = [];
   opt = merge_options(opt, varargin{:});
   
   num_timesteps = numel(schedule.step.val);
   num_inj_timesteps = numel(schedule.step.val(schedule.step.control==1));
   tSteps = opt.tStep;
   if isempty(tSteps)
      numSteps = numel(wellSols); % @@ use wellSols or schedule?
      tSteps = (1:numSteps)';
      dts = schedule.step.val;
   else
      assert(numel(tSteps) == 1);
      numSteps = 1;
      dts = schedule.step.val(opt.tStep);
   end
   
   obj = repmat({[]}, numSteps, 1);
   total_inj = zeros(numSteps, numel(wellSols{1,1}));
   obj_per_well = zeros(numSteps, numel(wellSols{1,1}));

   for step = 1:numSteps
      sol = wellSols{tSteps(step)};
      nW = numel(sol);
      pBHP = zeros(nW, 1); % place holder
      qGs = vertcat(sol.qGs);
      qWs = vertcat(sol.qWs);
      p = zeros(Gt.cells.num,1); % place holder
      sG = zeros(Gt.cells.num,1); % place holder
      if opt.ComputePartials
         %[p, sG, pBHP, qWs, qGs] = initVariablesADI(p, sG, pBHP, qWs, qGs);%#ok
         [p, sG, qWs, qGs, pBHP] = initVariablesADI(p, sG, qWs, qGs, pBHP);%#ok
      end
      dt = dts(step);
      injInx = (vertcat(sol.sign) > 0);
      
      
      % For each time step, compute the well cost
      
      %cost_for_having = ones(nW, 1) * (initial_cost/num_inj_timesteps); % USD
      cost_for_having = ones(nW, 1) * (initial_cost/num_inj_timesteps) .* max(0, sign(qGs - sqrt(eps))); % @@ assuming sqrt(eps) is min val
                    
      % NB: operation_cost is [USD/tonne] (1 tonne = 1e3 kg)
      % NB: dt*qGs is [m3] at ref depth
      % NB: dt*qGs*rhoGS is [kg]
      cost_for_operating = operation_rate * dt * injInx .* qGs .* fluid.rhoGS/1e3 .* max(0, sign(qGs - sqrt(eps)));     % USD
                       
      
      if ~nonlinear_well_cost
          
        % Linear total well cost:
        obj{step} = sum( (cost_for_having + cost_for_operating) .* max(0, sign(qGs - 2*sqrt(eps))) ); % USD
        % for plotting:
        obj_per_well(step,1:nW) = cost_for_having + cost_for_operating; % USD
      
      else
          
        % Non-linear total well cost:
        % Critical injection mass rate (for this time step) to make well's investment cost worthwhile:
        %q_crit = (initial_cost/num_inj_timesteps)/(operation_rate) * 1e3 / model.fluid.rhoGS / dt; % m3/s
        q_crit = (initial_cost/num_inj_timesteps)/(tax_credit_rate - operation_rate) * 1e3 / fluid.rhoGS / dt; % m3/s
        obj{step} = sum( cost_for_having .* tanh(qGs./(alpha*q_crit)) + cost_for_operating ); % USD
        % for plotting:
        obj_per_well(step,1:nW) = cost_for_having .* tanh(qGs./(alpha*q_crit)) + cost_for_operating; % USD
      
      end
      
      % for plotting purposes
      total_inj(step,1:nW) = qGs .* dt .* fluid.rhoGS .* max(0, sign(qGs - sqrt(eps))); % kg

      
      if (tSteps(step) == num_timesteps)
         if ~opt.ComputePartials
            fprintf('Well cost is %f USD. \n\n', sum(cell2mat(obj)) );
            %
            % total well cost
            if ~ishandle(12)
                figure(12)
            else
                % Avoid stealing focus if figure already exists
                set(0, 'CurrentFigure', 12); clf(12);
            end
            plot_total_well_cost_function( 12, initial_cost, operation_rate, tax_credit_rate, alpha ) % opens new figure
            hold on; plot(sum(total_inj)/1e9, sum(obj_per_well)/1e6, 'xk')
            text(sum(total_inj)/1e9, sum(obj_per_well)/1e6, cellstr(num2str([1:nW]')))
            %
            % max savings possible per well (assuming no leakage)
            if ~ishandle(13)
                figure(13)
            else
                % Avoid stealing focus if figure already exists
                set(0, 'CurrentFigure', 13); clf(13);
            end
            plot_total_savings_function( 13, initial_cost, operation_rate, tax_credit_rate, alpha )
            hold on; plot(sum(total_inj)/1e9, tax_credit_rate .* sum(total_inj)/1e9 - sum(obj_per_well)/1e6, 'xk')
            text(sum(total_inj)/1e9, tax_credit_rate .* sum(total_inj)/1e9 - sum(obj_per_well)/1e6, cellstr(num2str([1:nW]')))
         end
      end
   end

end


