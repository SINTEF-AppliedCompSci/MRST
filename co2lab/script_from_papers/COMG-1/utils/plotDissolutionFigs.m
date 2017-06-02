function plotDissolutionFigs(simres)

%% Produce panel plots for Figure 12

   produce_panel_plots(simres);


%% Produce detail plot for Figure 13

   produce_detail_plot(simres);

%% Produce plots for Figure 14

   produce_graphs(simres);

end

% ----------------------------------------------------------------------------

function produce_panel_plots(simres)

   fig_w = 1200;
   fig_h = 250;

   k     = 4; % case with residual saturation, dissolution and A = 0
   Gt    = simres{k}.Gt;
   fluid = setup_fluid(simres{k}.fluid_params);
   xc    = Gt.cells.centroids(:,1) / 1e3;
   zt    = zeros(size(xc));
   zb    = zt + Gt.cells.H;

   for nn = [40, 50, 70]
      f = figure(nn); clf;
      state = simres{k}.states{nn};
      sG_free = free_sg(state.s(:,2),state.sGmax, ...
                        struct('res_gas',fluid.res_gas, 'res_water', fluid.res_water));
      sGmax = state.sGmax;

      p = state.pressure; sG = state.s(:,2);
      h = (sG_free .* Gt.cells.H) ./ (1 - fluid.res_water);

      patch(xc([1 1:end end]), [zt(end) - 20; zt; zt(end) - 20],.7 * [1 1 1]);
      patch(xc([1 1:end end]), [zb(1) + 20; zb; zb(1) + 20],.7 * [1 1 1]);
      patch(xc([1:end end: - 1:1]),...
            [zt; zt(end: - 1:1)], myCOColor(1))
      patch(xc([1:end end: - 1:1]), [zt + h; zt(end: - 1:1)], myCOColor(2))

      h_max = (sGmax .* Gt.cells.H) ./ (1 - fluid.res_water);
      mm = minRs(p, sG, sGmax, fluid, Gt);
      h_res = Gt.cells.H .* ((1 - sG) .* state.rs - mm)/fluid.dis_max;


      assert(all(h_res > -1e-4))
      patch(xc([1:end end:-1:1]),...
            [zt + h; zt(end:-1:1) + h_max(end:-1:1)], myCOColor(3))
      patch(xc([1:end end:-1:1]),...
            [zt + h_max; zt(end:-1:1) + h_max(end:-1:1) + h_res(end:-1:1)], myCOColor(4))
      patch(xc([1:end end:-1:1]),...
            [zt + h_max + h_res; zb(end:-1:1)], myCOColor(5))

      line(xc([1:end]), zt + h_max, 'LineWidth', 2, 'Color', 'k');%#ok
      patch(xc([1 1:end end]), [zt(end) - 20; zt; zt(end) - 20], .7 * [1 1 1]);
      set(gca, 'YDir', 'reverse'), axis tight, axis([0 25 0 50]), axis off;
      pos = get(f, 'position');
      set(f, 'position', [pos(1), pos(2), fig_w, fig_h]);
   end
end

% ----------------------------------------------------------------------------

function produce_detail_plot(simres)

   k     = 8; % case with residual saturation, dissolution and A = 0
   Gt    = simres{8}.Gt;
   fluid = setup_fluid(simres{k}.fluid_params);
   xc    = Gt.cells.centroids(:,1) / 1e3;
   zt    = Gt.cells.z;
   zb    = zt + Gt.cells.H;

   z_tmp = max(Gt.cells.z) * ones(size(zt));
   for i = 2:numel(z_tmp)-1
      z_tmp(i) = max(zt(i:end));
   end
   z_tmp(end) = max(z_tmp(end-1), z_tmp(end));
   h_trap = z_tmp - Gt.cells.z;

   for nn=[50,80]
      figure(nn+100),clf,hold on;
      state = simres{k}.states{nn};
      sG_free = free_sg(state.s(:,2),state.sGmax, ...
                        struct('res_gas',fluid.res_gas, 'res_water', fluid.res_water));
      sGmax=state.sGmax;

      p = state.pressure;sG=state.s(:, 2);
      h = (sG_free .* Gt.cells.H) ./ (1 - fluid.res_water);

      patch(xc([1 1:end end]), [zt(end)-20; zt; zt(end)-20],.7*[1 1 1]);
      patch(xc([1 1:end end]), [zb(1)+20; zb; zb(1)+20],.7*[1 1 1]);
      patch(xc([1:end end:-1:1]), [zt + h_trap; zt(end:-1:1)], myCOColor(1))
      patch(xc([1:end end:-1:1]), [zt + h; zt(end:-1:1) + h_trap(end:-1:1)], myCOColor(2))

      h_max = (sGmax .* Gt.cells.H) ./ (1 - fluid.res_water);
      mm = minRs(p, sG, sGmax, fluid, Gt);
      h_res = Gt.cells.H .* ((1 - sG) .* state.rs - mm) / fluid.dis_max;

      assert(all(h_res>-1e-4))
      patch(xc([1:end end:-1:1]), [zt + h; zt(end:-1:1)+h_max(end:-1:1)],myCOColor(3))
      patch(xc([1:end end:-1:1]), ...
            [zt + h_max; zt(end:-1:1)+h_max(end:-1:1)+h_res(end:-1:1)], myCOColor(4))
      patch(xc([1:end end:-1:1]), [zt + h_max+h_res; zb(end:-1:1)], myCOColor(5))

      line(xc([1:end]), zt + h_max, 'LineWidth', 2, 'Color', 'k');%#ok
      patch(xc([1 1:end end]), [zt(end) - 20; zt; zt(end) - 20],.7 * [1 1 1]);
      set(gca, 'YDir', 'reverse'), axis tight

      % add innlet
      set(gca, 'FontSize', 19)
      set(gca, 'Ydir', 'reverse')
      xlabel('X km', 'FontSize', 16)
      ylabel('Depth m', 'FontSize', 16)
      set(gca, 'Color', 'none')
      h1 = gca;
      pos = get(h1, 'Position');%#ok
      h3 = copyobj(h1, gcf);

      set(h3, 'Position', get(h1, 'Position') * 0.3 + [0.55, 0.39, 0.07, 0])
      axes(h3);%#ok
      set(gca, 'FontSize', 16)
      set(gca, 'YDir', 'reverse');
      top = @(x) interp1(xc, zt, x);
      if(nn == 50)
         xx = [10, 10.5] + 2.5;
      elseif(nn == 80)
         xx = [21 23];
      else
         xx = [21 23];
      end
      matop = max(top(xx));
      mitop = min(top(xx));
      ax = [xx, mitop, matop];
      axis(ax);
      axes(h1);%#ok
      plot([ax(1), ax(2), ax(2), ax(1), ax(1)], [ax(3), ax(3), ax(4), ax(4), ax(3)], 'r')
      axes(h3);%#ok
      axis off

      xlabel(''), ylabel('')
      h4 = copyobj(h1, gcf);
      set(h4, 'Position', get(h1, 'Position') * 0.3 + [0.12, 0.635, 0.1, 0.0]);
      axes(h4)%#ok
      set(gca, 'YDir', 'reverse');
      if(nn == 50)
         xx = [5, 5.5];
      elseif(nn == 80)
         xx = [2.5 5];
      else
         xx = [2.5 5];
      end

      matop = max(top(xx));
      if(nn >= 80)
         matop = max(top(xx)) + 50;
      end
      mitop = min(top(xx));
      ax = [xx, mitop, matop];
      axis(ax);
      axes(h1);%#ok
      plot([ax(1), ax(2), ax(2), ax(1), ax(1)], [ax(3), ax(3), ax(4), ax(4), ax(3)], 'r')
      axes(h3)%#ok
      axes(h4)%#ok
      axis off
      %%

      h5 = copyobj(h1, gcf);
      set(h5, 'Position', get(h1, 'Position') * 0.3 + [0.30, 0.11, 0.3, 0.0]);
      axes(h5)%#ok
      set(gca, 'YDir', 'reverse');
      if(nn == 50)
         xx = [7, 8] + 2;
      elseif(nn == 80)
         xx = [7.03, 7.33] + 9.3;
      else
         xx = [7.03, 7.33] + 9.3;
      end

      matop = max(top(xx));
      if(nn == 80)
         matop = max(top(xx));
      end
      mitop = min(top(xx));
      ax = [xx, mitop, matop];
      axis(ax);
      axes(h1);%#ok
      plot([ax(1), ax(2), ax(2), ax(1), ax(1)], [ax(3), ax(3), ax(4), ax(4), ax(3)], 'r')
      axes(h3);%#ok
      axes(h4);%#ok
      axes(h5);%#ok
      axis off

      %%
      axx = {h3, h4, h5};
      for i = 1:numel(axx)
         axs = axx{i};
         ax = axis(axs);
         bpos = [ax([2, 4]), ax([2, 4]) - ax([1, 3])];
         xpos = mean(ax([1, 2]) - 0.03);
         bpos([1, 2]) = [xpos, top(xpos)];
         % pos5 = get(h5, 'Position');
         pos5 = dsxy2figxy_new(axs, bpos - [0 0 0 0]);
         pos1 = dsxy2figxy_new(h1, bpos + [0 0 0 0]);
         annotation('arrow', [pos1(1) pos5(1)], [pos1(2) pos5(2)], 'LineWidth',2)
      end
      drawnow;
   end
end

% ----------------------------------------------------------------------------

function produce_graphs(simres)

   legendtext = {'No dissolution (A=0)', 'Dissolution (A=0)', ...
                 'No dissolution (A=2)', 'Dissolution (A=2)'};
   linetype = {'b-', 'b--', 'r-', 'r--'};

   Gt = simres{4}.Gt;
   xc = Gt.cells.centroids(:, 1) / 1e3;

   for residual = [false, true]% residual saturation or not
      f = figure(); clf,
      set(f, 'Position', [0, 600, 700, 500]);
      kk = 1;
      for n = 1:2,% flat or non flat topsurface
         for dissolution = [false, true]
            k = 1 + (n - 1) * 4 + 2 * double(residual) + double(dissolution);
            state = simres{k}.states{70};
            fluid = setup_fluid(simres{k}.fluid_params);
            sG = free_sg(state.s(:, 2), state.sGmax,...
                         struct('res_gas', fluid.res_gas, 'res_water', fluid.res_water));
            hold on
            ff = filter_fun(n, xc);
            plot(xc, filter2(ff, sG .* Gt.columns.dz), linetype{kk}, 'LineWidth', 2);
            hold off
            drawnow;
            kk = kk + 1;
         end
      end
      axis tight
      ax = axis();
      axis([ 0 28 ax(3) ax(4)])
      set(gca, 'YDir', 'reverse', 'FontSize', 16);
      legend(legendtext{:}, 'Location', 'SouthWest');
   end
end

% ----------------------------------------------------------------------------
function ff = filter_fun(n, xc)
   if n == 1
      ff = 1;
   else
      xx = xc(150:650);
      ff = exp( -((xx - xc(400)) / (0.3)).^2);
      ff = ff / sum(ff);
   end
end

% ----------------------------------------------------------------------------

function fluid = setup_fluid(fluid_params)
   fluid = makeVEFluid(fluid_params.args{:});
   if fluid_params.is_instant
      fluid.dis_rate = 0;
   end
end
