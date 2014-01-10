function [state, report] = run_impes(fn)
require eclipse blackoil impes

[G, rock, fluid, state, wells, Trans, TSTEP, deck] = initEclipseModel(fn);

%%
try
   [d, f] = fileparts(fn);
   smry = readEclipseSummaryUnFmt(fullfile(d, f));

   figure(2);
   clf,
   subplot(2,2,1);

   plot(smry.TIME(2:end), smry.WGOR(1,2:end),'--r');
   hold on;

   subplot(2,2,2);

   plot(smry.TIME(2:end), smry.WBHP(:,2:end),'--');
   hold on;

   TSTEP = convertFrom(diff(smry.TIME), day);
catch ME
   disp(ME.message);
end

%%
forces = {'wells', wells};
porvol = poreVolume(G, rock);

%--------------------------------------------------------------------------
%% Run model --------------------------------------------------------------
%
start  = deck.RUNSPEC.START;
nstep  = numel(TSTEP); T = 0;

fprintf('\n\n%s Start run   %3d/%3d %s (%s)\n\n\n', ...
      repmat('=', [1, 10]), 0, nstep,   ...
      repmat('=', [1, 10]), datestr(start, 29));
report = [];
dt     = inf;


for k = 1 : nstep,

   fprintf('\n\n%s Report Step %3d/%3d %s (%s, %d days)\n\n\n', ...
           repmat('=', [1, 10]), k, nstep,   ...
           repmat('=', [1, 10]), datestr(start+sum(TSTEP(1:k))/day, 29),...
           TSTEP(k)/day);


   % reporting time step read from deck
   DT = TSTEP(k);

   % use reporting time step or previously adjusted (mini) time step
   dt = min(DT, dt);

   % t measures progress in current reporting time step [0, DT]
   t  = 0;

   done = false;
   while (~done)

      % flag to indicate if time step was acceptable
      ok = false;

      try

         fprintf('MINISTEP size %3.2g days...\n', convertTo(min(dt, DT-t), day));
         state0 = state;

         % attempt to solve impes system and update masses
         [state1, ddt, report] = impesTPFA2                             ...
            (state0, G, Trans, fluid, min(dt, DT-t), porvol, forces{:},  ...
            'report', report, 'Atol', 1e-8, 'Rtol', 1e-11);

         % reservoir fluid volume
         [u0,u0,u0,u0] = fluid.pvt(state0.pressure, state0.z);
         [u,u,u,u] = fluid.pvt(state1.pressure, state1.z);

         dplim = [5.0,  20.0  50.0]*barsa;
         dslim = [0.05, 0.1, 0.15];

         % compute max change in saturation and in reservoir pressure
         ds = max(abs(u(:)-u0(:)));
         dp = max(abs(state1.pressure-state0.pressure));

         ok = dp < dplim(3) && ds < dslim(3);

         if ok,

            % time step was acceptable:

            t     = t + ddt;
            state = state1;

         else

            % time step failed, report cause:

            if ~(dp < dplim(3)),
               c = find(abs(state1.pressure-state0.pressure) > dplim(3), 1, 'first');
               message = sprintf(...
                'Pressure changed more than %3.1e (bar) in cell %d.\n', ...
                dplim(3), c);
            else
               du = abs(u-u0);
               c = find(any(bsxfun(@gt, du, dslim(3)), 2), 1, 'first');
               message = sprintf(...
                  'Saturation changed more than %2.3e in cell %d.\n', ...
                  dslim(3), c);
            end

            fprintf('Step failed: %s\n', message);
         end


         % increase or decrease time step by scaling
         if     ds < dslim(1), dts = dt * dslim(1)/ds;
         elseif ds > dslim(2), dts = dt * dslim(2)/ds;
         else                  dts = dt;
         end

         if      dp < dplim(1), dtp = dt * dplim(1)/dp;
         elseif  dp > dplim(2), dtp = dt * dplim(2)/dp;
         else                   dtp = dt;
         end

         dt = min(dts, dtp);

      catch ME

         % agressively reduce time step
         dt = dt/2;
         fprintf('Step failed: %s\n', ME.message);
         if strcmp(ME.stack(1).name, 'update_pressure'),
            % This is probably well bhp becoming negative...
            assert(false);
         end

      end

      % Finish reporting time step DT
      done = ok && (~(t<DT));


   end

   plotstep(state, fluid, G, report);

   % reporting time step finished, update time T
   T = T + DT;
end
end




function plotstep(x, fluid, G, report)
   figs = [1,2];
   h = ishandle(figs);

   for i = find(~h), figure(figs(i)), end

   set_figure = @(f) set(0, 'CurrentFigure', f);

   [u, u, u, u] = fluid.pvt(x.pressure, x.z);
   s = bsxfun(@rdivide, u, sum(u,2));

   components = {'water', 'oil', 'gas'};
   phases     = {'aqua', 'liquid', 'vapor'};
   if true,

      set_figure(1)
      clf

      subplot(4, 2, 1)

        plotCellData(G, x.pressure/barsa);
        title('Pressure');
        axis tight;
        colorbar

      for c = 1 : 3,
         subplot(4,2,1+2*c)

         plotCellData(G, x.z(:,c));
         axis tight;
         title(sprintf('%s surface volume', components{c}));
         colorbar
      end

      subplot(4,2,2)

        plotCellData(G, sum(u, 2) - 1);
        axis tight;
        title('Volume error');
        colorbar

      for c = 1 : 3,
         subplot(4,2,2+2*c)

         plotCellData(G, s(:,c));
         axis tight;
         title(sprintf('%s saturation', phases{c}));
         colorbar
      end

      drawnow
   end

   if true
      set_figure(2);%cla;

      time = convertTo([report.TIME], day);
      GOR  = convertTo([report.WGOR], 1000*ft^3/stb);
      BHP  = convertTo([report.WBHP], psia);
      GPR  = convertTo([report.WGPR], 1000*ft^3/day);
      OPR  = convertTo([report.WOPR], stb/day);

      subplot(2,2,1);
      plot(time, GOR);
      title('WGOR');

      subplot(2,2,2);
      plot(time, BHP);
      title('WBHP');

      subplot(2,2,3);
      plot(time, GPR);
      title('WGPR');

      subplot(2,2,4);
      plot(time, OPR);
      title('WOPR');

   end
end
