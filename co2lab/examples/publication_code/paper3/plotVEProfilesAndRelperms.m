% Check that path exists
if exist('figs') ~= 7%#ok
   mkdir('figs');
end


clear all;
param=[0.2,0.1];
for entrypress=[false,true]

   gravity on;

    src   = param(1,1); srw =param(1,2);
    T_ref = 277+1000*30/1e3;
    mu    = [6e-2*milli 8e-4]*Pascal*second;
    rho   = [760 1200] .* kilogram/meter^3;


    potl_hyst = (src>0);
    H = 25;
    sample = 100;
    G = cartGrid([sample 1 1], [100 100 H]);
    G = computeGeometry(G);
    Gt = topSurfaceGrid(G);
    rock = struct('perm', 1000 * milli * darcy * ones(G.cells.num, 1), ...
                  'poro', 0.4 * ones(G.cells.num, 1));
    rock.parent = rock;
    fluidADI = initSimpleADIFluid('mu', [mu(2) mu(2) mu(1)],...
                                  'rho', [rho(2) rho(2), rho(1)],...
                                  'n', [1 1 1]);
    wfields = {'krO' 'krW', 'krG', 'pcWG', 'pcOG'};
    for i = 1:numel(wfields)
       if(isfield(fluidADI, wfields{i}))
          fluidADI = rmfield(fluidADI,wfields{i});
        end
    end
    fluidADI.pvMultR = @(p) 1 + (1e-5 / barsa) * (p - 100 * barsa);
    fluidADI.bW = @(p) 1 + (4.3e-5 / barsa) * (p - 100 * barsa);
    fluidADI.BW = @(p) 1 ./ fluidADI.bW(p);
    fluidADI.bG = boCO2(T_ref, fluidADI.rhoGS);
    fluidADI.BG = @(p) 1 ./ fluidADI.bG(p);
    fluidADI.bG = @(p) 1 + 0 * p; fluidADI.BG = @(p)  1 + 0 * p;
    fluidADI.bW = @(p) 1 + 0 * p; fluidADI.BW = @(p)  1 + 0 * p;
    fluidADI.rhoGS = 600;fluidADI.rhoWS = 1000;

    %%
    drho    = 400;
    C       = 4 * max(Gt.cells.H) * 0.4 * drho * norm(gravity);
    alpha   = 0.5;
    beta    = 3;
    samples = 100;

    if(entrypress)
       table_co2_1d = makeVEtables('invPc3D', @(p) max((C ./ (p + C)).^(1 / alpha) - src, srw),...
                                   'is_kscaled' , false        ,...
                                   'kr3D'       , @(s) s.^beta ,...
                                   'drho'       , drho         ,...
                                   'Gt'         , Gt           ,...
                                   'samples'    , samples);
    else
       table_co2_1d = makeVEtables('invPc3D', @(p) max((C ./ (p + C)).^(1 / alpha) - src * 0.0, srw),...
                                   'is_kscaled' , false        ,...
                                   'kr3D'       , @(s) s.^beta ,...
                                   'drho'       , drho         ,...
                                   'Gt'         , Gt           ,...
                                   'samples'    , samples);
    end
    S_tab           = linspace(0, 1, 10)';
    S_tab_w         = S_tab;
    kr_tab_w        = S_tab_w;
    table_water_1d  = struct('S', 1 - S_tab_w, 'kr', kr_tab_w, 'h', []);
    fluid           = fluidADI;
    fluid.invPc3D   = table_co2_1d.invPc3D;
    fluid.kr3D      = table_co2_1d.kr3D;
    fluid.res_water = srw;
    fluid.res_gas   = src;
    fluid.name      = 'variable';
    fluid_org       = fluid;

    table_co2_1d = makeVEtables('invPc3D', @(p) max((C ./ (p + C)).^(1 / alpha) - src, srw),...
                                'is_kscaled' , false        ,...
                                'kr3D'       , @(s) s.^beta ,...
                                'drho'       , drho         ,...
                                'Gt'         , Gt           ,...
                                'samples'    ,samples);

    fluid           = fluidADI;
    fluid.invPc3D   = table_co2_1d.invPc3D;
    fluid.kr3D      = table_co2_1d.kr3D;
    fluid.res_water = srw;
    fluid.res_gas   = src;
    fluid.name      = 'constant';
    fluid_org2      = fluid;
    fluid           = {fluid_org, fluid_org2};


    %%
    Hr         = Gt.cells.H(1);
    hh         = linspace(0, Gt.cells.H(1) * 4, 100);
    [H, H_max] = meshgrid(hh, hh);
    H_max_in   = max(H, H_max);
    msize      = size(H);
    org_SH     = nan(msize);
    org_krH    = nan(msize);
    org_s_max  = nan(msize);
    sara_SH    = nan(msize);
    sara_krH   = nan(msize);
    sara_s_max = nan(msize);
    p_ref      = 200 * barsa;
    kscale     = false;

    for i = 1:msize(1)
       for j = 1:msize(2)
          [s, pc, kr, SH, krH, s_max, fval] = veRelpermTester(H(i, j), p_ref, fluid_org, Hr, ...
                                                            'samples'   , 1000          , ...
                                                            'hs_max'    , H_max_in(i,j) , ...
                                                            'kscale'    , kscale        , ...
                                                            'const_res' , false); %#ok
          org_SH(i, j)    = SH;
          org_krH(i, j)   = krH;
          org_s_max(i, j) = s_max;
          [s, pc, kr, SH, krH, s_max, fval] = veRelpermTester(H(i, j), p_ref, fluid_org2, Hr, ...
                                                            'samples'   , 1000           , ...
                                                            'hs_max'    , H_max_in(i, j) , ...
                                                            'kscale'    , kscale         , ...
                                                            'const_res' , true); %#ok
          sara_SH(i, j)    = SH;
          sara_krH(i, j)   = krH;
          sara_s_max(i, j) = s_max;
       end
    end
    %%
    for i = 1:msize(1)
       for j = i + 1:msize(2)
          org_s_max(i, j) = org_s_max(i, i);
          sara_s_max(i, j) = sara_s_max(i, i);
       end
    end
    H_ff = Gt.cells.H(1);
    figure(31)
    clf; mesh(org_SH / Hr, org_s_max, org_krH ./ H_ff); hold on; mesh(sara_SH / Hr, sara_s_max, sara_krH ./ H_ff); xlabel('S')
    figure(32)
    clf; mesh(org_SH / Hr, org_s_max, H); hold on; mesh(sara_SH / Hr, sara_s_max, H); xlabel('S')
    scatteredInterpolant = @(x, y, z, method, extrap) TriScatteredInterp(x, y, z, method);%#ok % ..., METHOD

    %%
    kr_org  = scatteredInterpolant(org_SH(:) / Hr, org_s_max(:), org_krH(:) ./ H_ff, 'linear', 'none');
    kr_sara = scatteredInterpolant(sara_SH(:) / Hr, sara_s_max(:), sara_krH(:) ./ H_ff, 'linear', 'none');
    pc_org  = scatteredInterpolant(org_SH(:) / Hr, org_s_max(:), H(:), 'linear', 'none');
    pc_sara = scatteredInterpolant(sara_SH(:) / Hr, sara_s_max(:), H(:), 'linear', 'none');

    %%
    SSmax_H_org      = scatteredInterpolant(org_SH(:) / Hr, org_s_max(:), H(:), 'linear', 'none');
    SSmax_H_max_org  = scatteredInterpolant(org_SH(:) / Hr, org_s_max(:), H_max(:), 'linear', 'none');
    SSmax_H_sara     = scatteredInterpolant(sara_SH(:) / Hr, sara_s_max(:), H(:), 'linear', 'none');
    SSmax_H_max_sara = scatteredInterpolant(sara_SH(:) / Hr,sara_s_max(:),H_max(:),'linear','none');

    %%
    if(entrypress)
       myname = 'entry';
    else
       myname = 'noentry';
    end
    figure(41), clf, hold on
    nn = 10000;
    s = linspace(0, 0.5, nn);
    ee = ones(size(s));
    for s_max = [0, 0.2, 0.4]
       if(s_max == 0)
          plot(s, kr_org(s, s_max * ee), s, kr_sara(s, s_max * ee), '-', 'Linewidth', 2 )
       else
          plot(s, kr_org(s, s_max * ee), s, kr_sara(s, s_max * ee), '--', 'Linewidth', 2 )
       end
    end
    set(gca, 'FontSize', 16)
    xlabel('s_{n}'), ylabel('kr_{n}', 'Rotation', 0)

    print(['figs/relperm_hyst',myname,'.eps'],'-depsc2')

    figure(42), clf, hold on
    nn = 10000;
    s = linspace(0, 0.5, nn);
    ee = ones(size(s));
    for s_max = [0, 0.2, 0.4]
       if(s_max == 0)
          plot(s, pc_org(s, s_max * ee), s, pc_sara(s, s_max * ee), '-', 'Linewidth', 2 )
       else
          plot(s, pc_org(s, s_max * ee), s, pc_sara(s, s_max * ee), '--', 'Linewidth', 2 )
       end
    end
    set(gca, 'FontSize', 16)
    xlabel('s_{n}'),ylabel('p_{c}/\Delta \rho g','Rotation',90)
    print(['figs/pc_hyst',myname,'.eps'],'-depsc2')

    %%
    H_res = Gt.cells.H(1);
    kscale = nan; % sqrt(rock.poro ./ (rock.perm)) * fluidADI.surface_tension;

    %%
    for kk = 1:2
       if(kk == 1)
          S = 0.2;
          S_max = 0.5;
          H_org = SSmax_H_org(S, S_max);
          H_org_max = SSmax_H_max_org(S, S_max);
          H_mod = SSmax_H_sara(S, S_max);
          H_mod_max = SSmax_H_max_sara(S, S_max);
          myname = ['S_', num2str(S * 10), '_Smax_', num2str(S_max * 10)];
       else
          H_org = 7; % 0.5 * sco_max * H_res;
          H_mod = 7; % 0.5 * sco_max * H_res;
          H_org_max = 15; % sco_max * H_res;
          H_mod_max = 15; % sco_max * H_res;
          myname = ['H_', num2str(floor(H_org)), '_Hmax_', num2str(floor(H_org_max))];
       end

       if(entrypress)
          myname = ['entry_', myname]; %#ok
       else
          myname = ['noentry_',myname];%#ok
       end
       [s,pc,kr,SH_org,krH, s_max_org,fval_org]= veRelpermTester(H_org, p_ref, fluid_org, H_res, ...
                                                         'samples'   , 1000      , ...
                                                         'hs_max'    , H_org_max , ...
                                                         'kscale'    , kscale    , ...
                                                         'const_res' , false); %#ok
       [s,pc,kr,SH_mod,krH, s_max_mod,fval_mod]= veRelpermTester(H_mod, p_ref, fluid_org2, H_res, ...
                                                         'samples'   , 1000      , ...
                                                         'hs_max'    , H_mod_max , ...
                                                         'kscale'    , kscale    , ...
                                                         'const_res' , true);
       if(kk == 1)
          assert(abs(S - SH_mod / H_res) < 1e-2);
          assert(abs(S - SH_org / H_res) < 1e-2);
          assert(abs(S_max - s_max_mod) < 1e-2);
          assert(abs(S_max - s_max_org) < 1e-2);
       end

       s_org = @(z) interp1(H_org - fval_org.h, 1 - fval_org.s_h, z, 'linear', 0);
       sres_org = @(z) interp1(H_org_max - fval_org.h_max, (1 - fval_org.s_hmax), z, 'linear', 0);
       s_mod_max = @(z) interp1(H_mod_max - fval_mod.h_max, (1 - fval_mod.s_hmax), z, 'linear', 0);
       s_mod = @(z) interp1(H_mod - fval_mod.h, 1 - fval_mod.s_h,z,'linear',0);

       %%
       z = linspace(0, 25, 1000);
       % plot our model
       figure(81), clf, hold on
       x1 = fluid{1}.res_gas * (sres_org(z) - s_org(z));
       y1 = z;
       patch([0, fluid{1}.res_gas * (sres_org(z) - s_org(z)), 0], [z(1), z, z(end)], myCOColor(3))
       x = fluid{1}.res_gas * (sres_org(z) - s_org(z)) + s_org(z);
       y = z;
       xx = [x, x1(end: - 1:1)];
       yy = [y, y1(end: - 1:1)];
       patch(xx, yy, myCOColor(2))
       xx = [x, 0, 1, 1];
       yy = [y, y(end), y(end), 0];
       patch(xx, yy, myCOColor(5))
       plot(sres_org(z), z, 'r-', 'LineWidth', 2)
       set(gca, 'FontSize', 16)
       ylabel('depth')
       xlabel('s_{n}')

       ah = H_org / 2;
       ofigx = [(sres_org(ah) - s_org(ah)) * fluid_org.res_gas, (sres_org(ah) - s_org(ah)) * fluid_org.res_gas + s_org(ah)];
       ofigy = [ ah, ah];
       [figx, figy] = dsxy2figxy_new(ofigx, H_res - ofigy);
       annotation('doublearrow', figx, figy)
       text(mean(ofigx) - 0.05, mean(ofigy) + 1.5,'s_{eff}','FontSize',16)

       ah = H_org + (H_org_max - H_org) / 3;
       ofigx = [0, sres_org(ah) * fluid_org.res_gas];
       ofigy = [ ah, ah];
       [figx, figy] = dsxy2figxy_new(ofigx, H_res - ofigy);
       annotation('doublearrow', figx, figy)
       text(mean(ofigx) - 0.03, mean(ofigy) + 1.5, 's_{n,r}','FontSize',16)

       ah = H_org + (H_org_max - H_org) / 8;
       ofigx = [0, sres_org(ah)];
       ofigy = [ ah, ah];
       [figx, figy] = dsxy2figxy_new(ofigx, H_res - ofigy);
       annotation('doublearrow', figx, figy)
       text(mean(ofigx) - 0.03, mean(ofigy) + 1.5, 's_{max}','FontSize',16)

       set(gca,'Ydir','reverse')

       %%
       print(['figs/orginal_reconstuct',myname,'.eps'],'-depsc2')

       % plot model used elsewere
       figure(82), clf, hold on
       x1 = fluid{1}.res_gas * (H_mod_max>z & z>H_mod);
       y1 = z;
       patch([0, x1, 0], [z(1), z, z(end)], myCOColor(3))
       %
       x = s_mod(z);
       y = z;
       xx = [0, x, 0]; % , x1(end: - 1:1)];
       yy = [0, y, 25]; % , y1(end: - 1:1)];
       patch(xx, yy, myCOColor(2))
       xx = [x + x1, 0, 1, 1];
       yy = [y, y(end), y(end), 0];
       patch(xx, yy, myCOColor(5))
       plot(s_mod_max(z), z, 'r-', 'LineWidth', 2)
       set(gca, 'Ydir', 'reverse')
       set(gca, 'FontSize', 16)
       ylabel('depth')
       xlabel('s_{n}')

       ah = H_mod / 2;
       ofigx = [0, s_mod(ah)];
       ofigy = [ ah, ah];
       [figx, figy] = dsxy2figxy_new(ofigx, H_res - ofigy);
       annotation('doublearrow', figx, figy)
       text(mean(ofigx) - 0.05, mean(ofigy) + 1.5, 's_{eff}', 'FontSize',16)

       ah = H_mod + (H_mod_max - H_mod) / 3;
       ofigx = [0, fluid_org.res_gas];
       ofigy = [ ah, ah];
       [figx, figy] = dsxy2figxy_new(ofigx, H_res - ofigy);
       annotation('doublearrow', figx, figy)
       text(mean(ofigx) - 0.03, mean(ofigy) + 1.5, 's_{n, r}','FontSize',16)

       ah = H_mod + (H_mod_max - H_mod) / 8;
       ofigx = [0, s_mod_max(ah)];
       ofigy = [ ah, ah];
       [figx, figy] = dsxy2figxy_new(ofigx, H_res - ofigy);
       annotation('doublearrow', figx, figy)
       text(mean(ofigx) - 0.03, mean(ofigy) + 1.5, 's_{max}','FontSize',16)


       print(['figs/mod_reconstuct', myname, '.eps'], '-depsc2')
    end
end
