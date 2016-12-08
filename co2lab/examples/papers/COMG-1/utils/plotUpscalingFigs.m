function plotUpscalingFigs(simres)

%% Produce left plot of figure 11 in paper
   legendtext = {'Fine scale', 'Accretion layer',...
                 'Analytic: square', 'Analytic: sinus'};

   linetype = {'k-', 'b-','r-','g-'};
   ff = exp(-linspace( -25, 25, 501).^2); ff = ff' / sum(ff);
   figure; hold on;
   for k = 1:4
      Gt = simres{k}.Gt;
      xc = Gt.cells.centroids(:, 1) / 1e3;
      state = simres{k}.states{end - 70};
      
      fluid_params = simres{k}.fluid_params;
      fluid = makeVEFluid(fluid_params.args{:});
      if fluid_params.is_instant
         fluid.dis_rate = 0;
      end

      sG = free_sg(state.s(:, 2), state.sGmax,...
                   struct('res_gas', fluid.res_gas, ...
                          'res_water', fluid.res_water));
      hold on
      plot(xc, filter2(ff, sG .* Gt.columns.dz), linetype{k}, 'LineWidth', 2);
   end
   axis tight
   set(gca, 'YDir', 'reverse', 'FontSize', 12);
   h = legend(legendtext{:}, 'Location', 'SouthEast'); set(h, 'FontSize', 14);
   set(gcf, 'Position', [680 580 800 420]);
end
