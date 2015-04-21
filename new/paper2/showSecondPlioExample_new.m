function showSecondPlioExample_new()

   moduleCheck('coarsegrid', 'deckformat', 'mex', 'ad-fi', 'ad-props');
   gravity reset on;
   
   %% Ascertain presence of saved results, and load them
   savefiles = {'data/secondPlioExample_1200_1.mat', ...
                'data/secondPlioExample_1200_2.mat', ...
                'data/secondPlioExample_1200_3.mat'};
   present = cellfun(@(file) exist(file, 'file'), savefiles);
   if (all(present) && ask_user('Saved result found. Re-use? [y/n] '))
      fprintf('Re-using old result.\n');
   else
      fprintf('Recomputing result. \n');
      secondPlioExample; % running script
   
   end   
   outcomes = cellfun(@(name) load(name), savefiles);
   Gt = outcomes(1).Gt;
   W  = outcomes(1).schedule.control(1).W;
   time = [0;cumsum(outcomes(1).schedule.step.val)/year];
   
   
   %% Pliocenesand: effect of dissolution
   
   leg_methods = {'Compressible', ...
                  'Instant dissolution','Dissolution rate', ''};
   models = [1,3,2];
   mysteps = [69 91];
   ptimes = [710 1510];
   leg_methods=leg_methods(models);
   m = numel(models);
   col = .5*(lines(m) + ones(m,3));
   x = Gt.cells.centroids(:,1);
   y = Gt.cells.centroids(:,2);
   for k = 1:numel(ptimes)
      [d, i] = min((time - ptimes(k)).^2); 
      mysteps(k) = i;
   end
   
   ptimes=time(mysteps);
   
   for step = mysteps
      figure; hold on;
      has_rs = []; 
      for n = 1:3
         time_tmp = [0; cumsum(res{models(n)}.control.step.val) / year]; 
         assert(time_tmp(step) == time(step));

         % calculate properties
         state=res{models(n)}.states{step};
         p= state.pressure;
         sG= state.s(:,2);
         sGmax=state.sGmax;
         sG = free_sg(sG, sGmax, ...
                      struct('res_gas',fluid.res_gas, 'res_water', ...
                             fluid.res_water));
         fluid = outcomes(models(m)).fluid;
      
         drho = fluid.rhoWS .* fluid.bW(p) - fluid.rhoGS .* fluid.bG(p); 
         h = (fluid.pcWG(sG, p, 'sGmax', sGmax)) ./ (drho * norm(gravity())); 
         h_max = (fluid.pcWG(sGmax, p, 'sGmax', sGmax)) ./ (drho * norm(gravity()));
         if(n == 3)
            mm = minRs(p, sG, sGmax, fluid, Gt); % .* Gt.cells.H; 
            h_res_diff = Gt.cells.H .* ((1 - sG) .* state.rs - mm) / fluid.dis_max; 
         else
            h_res_diff = 0 * Gt.cells.H;
         end

         assert(all(h_res_diff >= -1e-5)); 
         h_res = h_max + h_res_diff;  
      
         FF_cc = TriScatteredInterp(x, y, h_max); % #ok
         [xx, yy] = meshgrid(linspace(min(x), max(x), 200),...
                             linspace(min(y), max(y), 200)); 
         cc = FF_cc(xx, yy); 
         cc(cc<.1) = NaN; 
         [c, a] = contourf(xx, yy, cc, 1e-1); 
         set(get(a, 'Children'), 'FaceColor', col(n, :), 'EdgeColor', col(n, :)); 
         set(a, 'EdgeColor', col(n, :));
         
      end
   end
   
   %% Draw inventory plot

   leg_traps = {'Dissolved','Residual (traps)', 'Residual', ...
                'Residual (plume)', 'Movable (traps)', ...
                'Movable (plume)', 'Leaked'};
   
   
end


