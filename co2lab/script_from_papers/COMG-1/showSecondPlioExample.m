function showSecondPlioExample()

   mrstModule add coarsegrid deckformat mex ad-core ad-props
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
   outcomes = cellfun(@(name) load(name), savefiles, 'uniformoutput', false);
   Gt = outcomes{1}.Gt;
   W  = outcomes{1}.schedule.control(1).W;
   time = [0;cumsum(outcomes{1}.schedule.step.val)/year];
   rock2D = outcomes{1}.rock2D;


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
      [d, i] = min((time - ptimes(k)).^2);%#ok
      mysteps(k) = i;
   end

   for step = mysteps
      figure; hold on;

      for n = 1:3
         time_tmp = [0; cumsum(outcomes{models(n)}.schedule.step.val) / year];
         assert(time_tmp(step) == time(step));

         % calculate properties
         state=outcomes{models(n)}.states{step};
         p= state.pressure;
         sG= state.s(:,2);
         sGmax=state.sGmax;
         fluid = outcomes{models(n)}.fluid;

         sG = free_sg(sG, sGmax, ...
                      struct('res_gas',fluid.res_gas, 'res_water', ...
                             fluid.res_water));

         drho = fluid.rhoWS .* fluid.bW(p) - fluid.rhoGS .* fluid.bG(p);
         h = (fluid.pcWG(sG, p, 'sGmax', sGmax)) ./ (drho * norm(gravity()));%#ok
         h_max = (fluid.pcWG(sGmax, p, 'sGmax', sGmax)) ./ (drho * norm(gravity()));
         if(n == 3)
            mm = minRs(p, sG, sGmax, fluid, Gt); % .* Gt.cells.H;
            h_res_diff = Gt.cells.H .* ((1 - sG) .* state.rs - mm) / fluid.dis_max;
         else
            h_res_diff = 0 * Gt.cells.H;
         end

         assert(all(h_res_diff >= -1e-5));
         h_res = h_max + h_res_diff;%#ok

         FF_cc = TriScatteredInterp(x, y, h_max); %#ok
         [xx, yy] = meshgrid(linspace(min(x), max(x), 200),...
                             linspace(min(y), max(y), 200));
         cc = FF_cc(xx, yy);
         cc(cc<.1) = NaN;
         [c, a] = contourf(xx, yy, cc, 1e-1);%#ok
         if verLessThan('matlab', '8.4')
            set(get(a, 'Children'), 'FaceColor', col(n, :), 'EdgeColor', col(n, :));
            set(a, 'EdgeColor', col(n, :));
         else
            set(a, 'facecolor', col(n, :));
            set(a, 'edgecolor', col(n, :));
         end
      end
      FF_zz = TriScatteredInterp(x, y, Gt.cells.z); %#ok
      zvec = FF_zz(xx, yy);
      contour(xx, yy, zvec, 40, 'k');

      axis equal tight off;
      plot(Gt.cells.centroids(W.cells, 1),...
           Gt.cells.centroids(W.cells, 2), 'ok',...
           'MarkerSize', 10, 'MarkerFaceColor', 'black')
      [legh, objh, outh, outm] = legend(leg_methods, 'location', 'east');%#ok
      set(legh, 'Fontsize', 12);
      for n = m + 1:2 * m
         set(get(get(objh(n), 'Children'), 'Children'), 'LineWidth', 4);
      end
      title(['after ' num2str(floor(time(step))), ' years'], 'FontSize',12)
   end

   %% Draw inventory plot

   leg_traps = {'Dissolved','Residual (traps)', 'Residual', ...
                'Residual (plume)', 'Movable (traps)', ...
                'Movable (plume)', 'Leaked'};
   ta = trapAnalysis(Gt, false);

   % Calculating all masses
   for m=1:3
      nsteps = numel(outcomes{m}.states);
      outcomes{m}.tot_masses = nan(nsteps, 3); % #ok
      outcomes{m}.masses = nan(nsteps, 8); % #ok
      schedule = outcomes{m}.schedule;
      for step=1:numel(outcomes{m}.states)
         % calculate properties
         state=outcomes{m}.states{step};
         p= state.pressure;
         sG= state.s(:,2);
         sGmax=state.sGmax;
         sG_free = free_sg(sG, sGmax, ...
                           struct('res_gas',fluid.res_gas, 'res_water', fluid.res_water));
         if(m == 2)
            sGmax = -((1 - fluid.res_water - fluid.res_gas) * sG_free - ...
                      (1 - fluid.res_water) * sG) ./ fluid.res_gas;
        end
        assert(all(sGmax>=0));

        fluid = outcomes{m}.fluid;
        drho = fluid.rhoWS .* fluid.bW(p) - fluid.rhoGS .* fluid.bG(p);
        h = (fluid.pcWG(sG, p, 'sGmax', sGmax)) ./ (drho * norm(gravity()));
        h_max = (fluid.pcWG(sGmax, p, 'sGmax', sGmax)) ./ (drho *norm(gravity()));
        if(n == 3)
           mm = minRs(p, sG, sGmax, fluid, Gt); % .* Gt.cells.H;
           h_res_diff = Gt.cells.H .* ((1 - sG) .* state.rs - mm) / fluid.dis_max;
        else
           h_res_diff = 0 * Gt.cells.H;
        end
        h_max_tmp = (sGmax / (1 - fluid.res_water)) .* Gt.cells.H;
        h_tmp = (sG_free / (1 - fluid.res_water)) .* Gt.cells.H;

        ss = (h * (1 - fluid.res_water) + (h_max - h) * fluid.res_gas) ./ Gt.cells.H;

        assert(all(abs(state.s(:, 2) - ss)<1e-3))
        assert(all(abs(h_max - h_max_tmp)<1e-3)); % valid for sharp interface when mass calculation is valid
        assert(all(abs(h - h_tmp)<1e-3)); % valid for sharp interface when mass calculation is valid

        assert(all(h_res_diff >= -1e-5));
        h_res = h_max + h_res_diff;%#ok
        % calulate masses for later
        if(step == 1)
           totMass = 0;
        else
           if(schedule.step.control(step - 1) == 1)
              dT = schedule.step.val(step - 1);
              totMass = totMass + fluid.rhoGS * W(1).val * dT;
           else
              % is no well
           end
        end

        masses = phaseMassesVEADI(Gt, state, rock2D, fluid);
        outcomes{m}.tot_masses(step, :) = masses; % #ok
        co2mass = masses(1) + masses(3);%#ok
        % h and h_max;
        state.h = h;
        state.h_max = h_max;
        % calculated distributions only valid for sharp interface.
        dh = Gt.cells.z * 0; % no subscale trapping
        if isfield(state, 'rs')
           rs = state.rs;
        else
           rs = 0;
        end
        masses = massTrappingDistributionVEADI(Gt, state.pressure, state.s(:,2), ...
                                               state.s(:,1), h, h_max, rock2D, ...
                                               fluid, ta, dh, 'rs', rs);
        
        % store all masses as a matrix with rows repersenting time
        outcomes{m}.masses(step, :) = [masses, totMass]; % #ok
      end
   end

% Draw inventory plot

figure(); set(gcf, 'Position', [100 100 1400 600]);
xw = 0.8 / numel(models);
xa = 0.2 / (numel(models) + 1);
for i = 1:numel(models)
   ind = [1:5, 7, 8];
   masses = outcomes{models(i)}.masses(:, ind);
   masses(:, end) = masses(:, end) - sum(masses(:, 1:end - 1), 2);
   axes('position', [(i - 1) * xw + i * xa,.17, xw,.78]); %#ok
   set(gca, 'FontSize', 14);
   hp = area(time, masses / 1e9);
   col = getInventoryColors(1:7);
   for j = 1:numel(hp), set(hp(j), 'FaceColor', col(j, :)); end
   axis tight
   title(leg_methods{i});

end
axes('position', [0.05 0.01 0.9 0.03], 'visible', 'off');
hL = legend(hp, leg_traps);
set(hL, 'orientation', 'horizontal', 'position',[0.05 0.02 0.9 0.03]);

end
