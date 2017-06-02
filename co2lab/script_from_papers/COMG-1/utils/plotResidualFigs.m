function plotResidualFigs(simres)

%% Produce plots for figure 5 and 6 in paper

   legendtext = {'No residual (A=0)', 'Residual (A=0)', ...
                 'No residual (A=2)', 'Residual (A=2)'};

   xc = simres{1}.Gt.cells.centroids(:, 1) / 1e3;
   ffun = @(num) filter_fun(num, xc);

   %% First figure
   figure; hold on;

   ff = ffun(1);
   for k = 1:2
      fluid_params = simres{k}.fluid_params;
      fluid = makeVEFluid(fluid_params.args{:});
      if fluid_params.is_instant
         fluid.dis_rate = 0;
      end
      
      state = simres{k}.states{end-70};
      linetype = {'b-', 'r-'};
      sG = free_sg(state.s(:, 2), state.sGmax ,...
                   struct('res_gas', fluid.res_gas, ...
                          'res_water', fluid.res_water));
      plot(xc, filter2(ff, sG .* simres{k}.Gt.columns.dz), linetype{k}, 'LineWidth', 2);
   end
   axis tight
   set(gca, 'YDir', 'reverse', 'FontSize', 16);
   legend(legendtext{1:2}, 'Location', 'SouthEast');

   %% Second figure
   figure; hold on;

   for k = 1:4
      fluid_params = simres{k}.fluid_params;
      fluid = makeVEFluid(fluid_params.args{:});
      if fluid_params.is_instant
         fluid.dis_rate = 0;
      end
      
      state = simres{k}.states{end-70};
      linetype = {'b--', 'r--', 'b-', 'r-'};
      ff = ffun(k);
      sG = free_sg(state.s(:, 2), state.sGmax ,...
                   struct('res_gas', fluid.res_gas, ...
                          'res_water', fluid.res_water));
      plot(xc, filter2(ff, sG .* simres{k}.Gt.columns.dz), linetype{k}, 'LineWidth', 2);
   end
   axis tight
   set(gca, 'YDir', 'reverse', 'FontSize', 16);
   legend(legendtext{:}, 'location','east');

   %% Third figure (solution in physical space)

   figure
   set(gcf,'PaperPositionMode','auto');
   p = get(gcf, 'Position'); set(gcf, 'Position', [p(1:2) 900 300]);

   % Main plot
   state = simres{3}.states{end - 70};
   Gt = simres{3}.Gt;
   ff = ffun(3);
   sG = free_sg(state.s(:, 2), state.sGmax, struct('res_gas', 0, 'res_water', 0));
   z1 = sG .* Gt.columns.dz;
   z2 = filter2(ff, z1);
   hold on
   plot(xc, z1, 'k-');
   plot(xc, z2, 'k-', 'LineWidth', 1);
   plot([15 16 16 15 15]', [0.5 0.5 5 5 0.5], 'r', 'LineWidth', 1);
   hold off
   h1 = gca;
   set(gca, 'YDir', 'reverse');

   % Inlet: zoom of subscale solution and average
   i = (xc >= 15) & (xc <= 16);
   axes('Position', get(h1, 'Position') * 0.3 + [0.6, 0.6, 0, 0]);
   plot(xc(i), z1(i), 'k-');
   hold on; plot(xc(i), z2(i), 'k-', 'LineWidth', 2); hold off
   set(gca, 'YDir', 'reverse'); axis tight


   % Inlet: zoom of the fluid distribution
   x = xc(i);
   zt = Gt.cells.z(i);
   zc = zt + z1(i);
   h4 = axes('Position', get(h1, 'Position') * 0.3 + [0.6, 0.2, 0, 0]);%#ok
   patch([x; x(end:-1:1)], [zt; zt(end:-1:1) + 50], myCOColor(5));
   patch([x; x(end:-1:1)], [zt; zc(end:-1:1)], myCOColor(2));
   set(gca, 'YDir', 'reverse'); axis tight
   set(gca, 'YLim', min(zt) + [0 25]);

end

% ----------------------------------------------------------------------------

function ff = filter_fun(k, xc)
   if k < 3
      ff = 1;
   else
      xx = xc(150:650);
      ff = exp( -((xx - xc(400)) / (0.3)).^2);
      ff = ff / sum(ff);
   end
end
