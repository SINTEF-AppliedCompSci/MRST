function plotTopSurfaceExample(simres)

   gravity on;

%% Produce plots for figure 11

   % Simulation indices for simulations with dissolution, caprock
   % oscillations, residual saturation and one of the fluid types 'sharp
   % interface', 'linear cap.' and 'P-scaled table'.
   ix = [10 12 14];
   top_surface_reconstruction(simres(ix));

%% Produce graphs for figure 12

   plot_plume_shapes(simres);

end

% ----------------------------------------------------------------------------

function plot_plume_shapes(res)

   legendtext = {'No dissolution (A=0)', 'Dissolution (A=0)', ...
                 'No dissolution (A=2)', 'Dissolution (A=2)'};

   line_colors={'k','b','m','g'};

   Gt_rough = res{9}.Gt;
   Gt       = res{1}.Gt;

   z = Gt_rough.cells.z;
   zt = max(z) * ones(size(z));
   for i = 2:numel(zt) - 1
      zt(i) = max(z(i:end));
   end
   zt(end) = max(zt(end - 1), z(end));
   ht = zt - z;
   ff = exp( - linspace( - 25, 25, 501).^2); ff = ff' / sum(ff);
   hts = filter2(ff, ht);

   % Now setting z to the smooth version of the caprock

   fluid_types = {'sharp interface','linear cap','P-scaled table', ...
                  'P-K-scaled table'};
   count = 1;
   for smooth = [true false]

      f = figure(); clf; hold on;
      set(f, 'Position', [0, 600, 700, 500]);

      for i = 1:4 % fluid type

         for dissolution = [false, true]

            state = res{count}.states{end-140};
            fluid = setup_fluid(res{count}.fluid_params);
            sG = free_sg(state.s(:,2),state.sGmax, ...
                         struct('res_gas',fluid.res_gas, ...
                                'res_water', fluid.res_water));
            hold on;
            if (dissolution)
               mline = [line_colors{i}, '--'];
            else
               mline = [line_colors{i}, '-'];
            end

            xc = Gt.cells.centroids(:,1)/1e3;
            xx = xc(150:650);
            ff = exp(-((xx-xc(400))/(0.3)).^2);
            ff = ff/sum(ff);

            plot(xc, filter2(ff, sG .* Gt.cells.H), mline, 'LineWidth', 2);

            drawnow;
            count = count + 1;
         end

      end
      if ~smooth
         plot(xc(50:end-50), hts(50:end-50), 'r-', 'LineWidth', 2);
      end

      axis tight;
      ax = axis();%#ok
      set(gca, 'YDir', 'reverse', 'FontSize', 16);
      % Add to legends one for color and one for linestyle
      ivec = [1, 2, 4, 3];
      if smooth
         addLegends(gca, fluid_types(ivec), line_colors(ivec), legendtext(1:2), ...
                    {' - ', ' -- '});
      else
         addLegends(gca, fluid_types(ivec), line_colors(ivec), legendtext(3:4), ...
                    {' - ', ' -- '});
      end
   end
end

% ----------------------------------------------------------------------------

function top_surface_reconstruction(res)

   Gt = res{1}.Gt;
   z  = Gt.cells.z;
   zt = max(z)*ones(size(z));
   for i=2:numel(zt)-1
      zt(i)=max(z(i:end));
   end
   zt(end) = max(zt(end-1),z(end));
   ht  = zt - z;
   h_trap = ht;

   rock = struct('perm', 100 * milli * darcy * ones(Gt.parent.cells.num, 1),...
                 'poro', 0.2 * ones(Gt.parent.cells.num, 1));
   rock2D = averageRock(rock, Gt);



   for f_ix = 1:3 % fluid type

      fluid = setup_fluid(res{f_ix}.fluid_params);
      Gt    = res{f_ix}.Gt;

      for nn = [100 160] % time step
         state = res{f_ix}.states{nn};

         sG_free = free_sg(state.s(:,2),state.sGmax, ...
                           struct('res_gas',fluid.res_gas, ...
                                  'res_water', fluid.res_water));

         f = figure(); clf, hold on;
         set(f, 'Position', [0, 600, 700, 500]);
         p          = state.pressure;
         sG         = state.s(:,2);
         sGmax      = state.sGmax;
         drho       = fluid.rhoWS .* fluid.bW(p) - fluid.rhoGS .* fluid.bG(p);
         h          = (fluid.pcWG(sG, p, 'sGmax', sGmax)) ./ (drho * norm(gravity()));
         h_max      = (fluid.pcWG(sGmax, p, 'sGmax', sGmax)) ./ (drho * norm(gravity()));
         mm         = minRs(p, sG, sGmax, fluid, Gt);
         h_res_diff = Gt.cells.H .* ((1 - sG) .* state.rs - mm) / fluid.dis_max;

         assert(all(h_res_diff >= -1e-5));
         h_res = h_max + h_res_diff;

         zt = Gt.cells.z;
         zb = zt + Gt.cells.H;
         xc = Gt.cells.centroids(:,1)/1e3;

         patch(xc([1 1:end end]), [zt(end) - 20; zt; zt(end) - 20],.7 * [1 1 1]);
         patch(xc([1 1:end end]), [zb(1) + 20; zb; zb(1) + 20],.7 * [1 1 1]);

         patch(xc([1:end end: - 1:1]),...
               [zt + h_trap; zt(end: - 1:1)], myCOColor(1))
         patch(xc([1:end end: - 1:1]),...
               [zt + h; zt(end: - 1:1) + h_trap(end: - 1:1)], myCOColor(2))
         patch(xc([1:end end: - 1:1]),...
               [zt + h; zt(end: - 1:1) + h_max(end: - 1:1)], myCOColor(3))
         patch(xc([1:end end: - 1:1]),...
               [zt + h_max; zt(end: - 1:1) + h_res(end: - 1:1)], myCOColor(4))
         patch(xc([1:end end: - 1:1]),...
               [zt + h_res; zb(end: - 1:1)], myCOColor(5))
         assert(all(h_max >= h));
         % line(xc([1:end]), zt + h_max, 'LineWidth', 2, 'Color', 'k')
         patch(xc([1 1:end end]), [zt(end) - 20; zt; zt(end) - 20], 0.7 * [1 1 1]);

         set(gca, 'YDir', 'reverse'), axis tight

         drawnow;
         hold on

         set(gca, 'FontSize', 19)
         set(gca, 'Ydir', 'reverse')
         xlabel('X km', 'FontSize', 16)
         ylabel('Depth m', 'FontSize', 16)
         set(gca, 'Color', 'none')

         h1 = gca;
         if(nn == 50)
            cells = [89, 135, 200, 270, 89];
         else% (nn == 80)
            cells = [89, 270, 350, 405, 89];
         end

         for ii = 1:numel(cells)
            cell = cells(ii);
            ha{ii} = axes();%#ok
            if(ii == 1)
               set(ha{ii}, 'Position', get(h1, 'Position') * 0.3 + [0.4, 0.78,- 0.05,- 0.15])
            elseif(ii == 2)
               set(ha{ii}, 'Position', get(h1, 'Position') * 0.3 + [0.30, 0.13, 0.0, 0.0]);
            elseif(ii == 3)
               set(ha{ii}, 'Position', get(h1, 'Position') * 0.3 + [0.60, 0.13, 0.0, 0.0]);
            elseif(ii == 4)
               set(ha{ii}, 'Position', get(h1, 'Position') * 0.3 + [0.60, 0.44, 0.00, 0])
            elseif(ii == 5)
               set(ha{ii}, 'Position', get(h1, 'Position') * 0.3 + [0.15, 0.635, 0.0, 0.0]);
            else
               error('Wrong cell');
            end

            hold on
            rock = rock2D;
            kscale = sqrt(rock.poro(cell) ./ (rock.perm(cell))) * fluid.surface_tension;
            H = Gt.cells.H(cell);

            h_int = h(cell);
            h_int_max = h_max(cell);
            h_int_res_b = h_res(cell);

            ss = []; hh = []; krk = []; ss_max = []; hh_max = [];

            [s, pc, kr, s_max, s_free, fval] = veRelpermTesterCell(h_int, drho(cell), fluid, H, ....
                                                              'samples', 1000, ...
                                                              'hs_max', h_int_max, ...
                                                              'kscale', kscale);
            mtol = 1e-2;
            assert(abs(s_free - sG_free(cell))<mtol);
            assert(abs(s - sG(cell))<mtol);
            assert(abs(s_max - sGmax(cell))<mtol);
            pc_cells = fluid.pcWG(sG, p, 'sGmax', sGmax);
            assert(abs(pc - pc_cells(cell))<mtol);
            krG = fluid.krG(sG, p, 'sGmax', sGmax);
            assert(abs(kr - krG(cell))<mtol);

            ss = [ss, 1 - fval.s_h];%#ok
            ss_max = [ss_max, 1 - fval.s_hmax];%#ok
            hh_max = [hh_max, h_int_max - fval.h_max];%#ok
            krk = [krk, fval.kr_h];%#ok
            hh = [hh, h_int - fval.h];%#ok

            s_max = @(d) interp1(hh_max, ss_max, d, 'linear', 0);
            if(h_int>0)
               s = @(d) interp1(hh, ss, d, 'linear',0);
            else
               s = @(d) d*0;
            end

            res_s = @(d) double(d<(h_int_res_b));
            w_s = @(d) ones(size(d));
            d = linspace(0, 50, 1345)';

            ptch =@(x,y) [x;y(end:-1:1);x(1)];
            axes(ha{ii});%#ok

            assert(all(s(d) >= 0))
            s_res = (s_max(d) - s(d)) .* fluid.res_gas; % 'non flowing co2'
            patch(ptch(s_res, s(d) * 0), ptch(d, d), myCOColor(3))
            s_real = s(d) + s_res; % total saturation
            ind = find(d>h_trap(cell), 1, 'first');
            patch(ptch(s_real(1:ind), s_res(1:ind)), ptch(d(1:ind), d(1:ind)), myCOColor(1))
            patch(ptch(s_real(ind + 1:end), s_res(ind + 1:end)), ptch(d(ind + 1:end), d(ind + 1:end)), myCOColor(2))

            fix = res_s(d);
            patch(ptch(fix, s_real), ptch(d, d), myCOColor(4))
            patch(ptch(w_s(d), fix),ptch(d,d),myCOColor(5))

            line([0.0 1.0],[h_trap(cell) h_trap(cell)],'Color',myCOColor(1),'LineWidth',1);

            ss_max = s_max(d);
            ss_max(d>h_int_max) = nan;
            ind = find(d>h_int_max, 1, 'first');
            ss_max(ind) = 0;
            plot(ss_max, d, 'k'); % , 'LineWidth', 2)

            set(gca, 'YDir', 'reverse')
            box on
            set(gca, 'LineWidth', 2)
            set(gca, 'FontSize', 16)
            axes(h1);%#ok
            line([xc(cell), xc(cell)], [zt(cell) - 50, zb(cell) + 50], 'LineWidth', 4, 'Color', 'r')
            axes(ha{ii})%#ok
         end
         top = @(x) interp1(xc, zt, x);
         for i = [5 1 2 3 4]
            axs = ha{i};
            if(i<5)
               xpos = xc(cells(i));
               bpos = [xpos, top(xpos) + 25, 1 1];
               pos5 = dsxy2figxy_new(axs, [0 25 1 1]);
               pos1 = dsxy2figxy_new(h1, bpos + [0 0 0 0]);
               annotation('arrow', [pos1(1) pos5(1)], [pos1(2) pos5(2)], 'LineWidth', 2)
            else
               axes(ha{i});%#ok
               ax = axis();%#ok
               axis([0 1, 0, h_trap(cells(i)) * 5]);
            end
            axes(ha{i});%#ok
            axis off;
         end
      end
   end
end

function fluid = setup_fluid(fluid_params)
   fluid = makeVEFluid(fluid_params.args{:});
   if fluid_params.is_instant
      fluid.dis_rate = 0;
   end
end
