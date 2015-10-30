function exploreHammerfest( varargin )
% Explore CO2 storage capacity in Hammerfest Basin Aquifer.
%
% SYNOPSIS:
%   exploreHammerfest();
%   exploreHammerfest('coarsening',N);
%
%
% DESCRIPTION:
%   The Hammerfest Basin aquifer is made up of three layered formations:
%   Sto, Nordmela, and Tubaen. Nordmela is treated as a sealing formation
%   due to its low permeability, thus its storage capacity is not
%   evaluated.
%
%   Storage potential of Sto and Tubaen are assessed using trapping
%   analysis and VE simulation. Structural traps identified using
%   trapAnalysis are compared to the "Greater Structural Closures" and
%   "Prospects" identified by NPD in the CO2 Atlas.
%
%   One of the structural closures identified is the "Greater Snohvit
%   Field", the location of the Snohvit CO2 injection project.
%   Reservoir-scale simulation of Snohvit is impractical due to the low
%   resolution of the Hammerfest formations. Thus, we only assess the
%   storage potential of Hammerfest at the basin-scale.
%
%
% PARAMETERS:
%   'pn'/pv - List of optional property names/property values:
%                   
%    - coarsening: Coarsening factor. If set to one, a grid with
%             approximately one cell per datapoint is produced. If set to
%             two, every second datapoint in x and y direction is used,
%             giving a reduction to 1/4th size. Default: 1
%
% RETURNS:
%
%
% SEE ALSO:
%   exploreCapacity, exploreSimulation


opt.coarsening = 4;
opt = merge_options(opt, varargin{:});


%% Construct top grids and rock properties

N = opt.coarsening;
[Gt_st, rock2D_st, petrodata_st] = getFormationTopGrid('Stofm',N);
[Gt_nd, rock2D_nd, petrodata_nd] = getFormationTopGrid('Nordmelafm',N);
[Gt_tu, rock2D_tu, petrodata_tu] = getFormationTopGrid('Tubaenfm',N);

% Plots
plotHammerfestLayers3D(Gt_st, Gt_nd, Gt_tu);
plotHammerfestRockProps(Gt_st, rock2D_st, Gt_nd, rock2D_nd, Gt_tu, rock2D_tu);

% save figure (optional)
%export_fig(gcf, [figDirName '/' 'FmRockProperties'], '-png','-transparent')


%% Perform trapping analysis
% NB: interactiveTrapping() can be called apart from exploreHammerfest(),
% however interactiveTrapping() does not access any heterogeneous rock
% properties
ta_st = trapAnalysis(Gt_st, false);
ta_tu = trapAnalysis(Gt_tu, false);




%% Injection Set-up and Simulation:

% Run injection scenario using entire Sto formation
[wellSols_st, states_st, sim_report_st, opt_st, var_st ] = ...
    runSnohvitInjectionScenario( Gt_st, rock2D_st );


% Alternatively, various injection locations can be specified here and
% passed into runSnohvitInjectionScenario()

% Physical coordinate(s):
wellCoords = [9.225e5, 7.988e6; 9.225e5 + 300, 7.988e6 + 300; ...
    9.225e5 + 400, 7.988e6 + 400; 9.225e5 + 500, 7.988e6 + 500];



% Run injection scenario using entire Sto formation
[wellSols_st, states_st, sim_report_st, opt_st, var_st ] = ...
    runSnohvitInjectionScenario( Gt_st, rock2D_st, 'wellCoords',wellCoords);


% Plot
plotCO2footprint(Gt_st, ta_st, states_st, opt_st, var_st);


% Analyze simulated injection/migration:
figure;
plotToolbar(Gt_st, states)

%%  
% Note that injection wells could be placed in other locations of the Sto
% formation, to exploit the full capacity of the structural traps. Also, it
% is quite clear by the outline of some structural traps where some of the
% NPD's prospects are location. Trapping analysis should find these
% prospect areas and display detail such as rock volume, pore volume,
% storage capacity, etc.

% could also use:
exploreCapacity( 'default_formation',  'Stofm',     ...
                 'grid_coarsening',     1,          ...
                 'seafloor_depth',      330*meter,  ...
                 'seafloor_temp',       4,          ...
                 'temp_gradient',       40           );
             
exploreSimulation(  'default_formation',  'Tubaenfm',     ...
                    'grid_coarsening',     3,          ...
                    'seafloor_depth',      330*meter,  ...
                    'seafloor_temp',       4,          ...
                    'temp_gradient',       40,         ...
                    'inj_time',            30 * year,  ...
                    'inj_steps',           30,         ...
                    'mig_time',            0 * year,   ...
                    'mig_steps',           0,          ...
                    'savefile',            'testsave'   );


% -------------------------------------------------------------------------
% ----------------------- DEPENDENT HELPER FUNCTIONS    -------------------
% -------------------------------------------------------------------------


end

% -------------------------------------------------------------------------
% -----------------------      HELPER FUNCTIONS    ------------------------
% -------------------------------------------------------------------------

function cellIndex = getCellIndex(Gt, Xcoord, Ycoord)
% Closest cell index of grid Gt corresponding to physical coordinate (X,Y)

    dv        = bsxfun(@minus, Gt.cells.centroids(:,1:2), [Xcoord, Ycoord]);
    [v, ind]  = min(sum(dv.^2, 2));
    cellIndex = ind; 
    % or Gt.cells.indexMap(i);
    
end
    
% -------------------------------------------------------------------------

function plotCO2footprint(Gt, ta, states, var, opt)
% Superimpose CO2 mass on contours, traps, etc. plotted by mapPlot

    figure;
    hold on
    subplot(2,2,[1 3])

    mapPlot(gcf, Gt, 'traps', ta.traps, 'trapalpha', 0.2, ...
        'rivers', ta.cell_lines, 'rivercolor', [1 0 0], ...
        'maplines', 20, 'wellcells',var.wellCellIndex);

    subplot(2,2,[2 4])

    CO2sat = states{end}.s(:,2);
    plumeOutline_SatTol = (0.01/100); % adjust this value if patch error occurs

    % to ensure grid data plotted above mapPlot contours
    Gt_tmp = Gt; Gt_tmp.nodes.z = -100*ones(Gt_tmp.nodes.num,1);
    plotCellData(Gt_tmp, CO2sat, CO2sat>plumeOutline_SatTol, 'EdgeColor','none')

    mapPlot(gcf, Gt, 'traps', ta.traps, 'trapalpha', 0.2, ...
        'rivers', ta.cell_lines, 'rivercolor', [1 0 0], ...
        'maplines', 20, 'wellcells',var.wellCellIndex);
    % caution: adding this map changes values of colorbar, thus, it is
    % important to adjust the colorbar afterwards.

    % We then re-plot the contours of the grid topology, using contour3()
    % to plot the contours at the elevation determined by function inside
    % drawContours3(). To ensure these final contours are on top of all
    % other plots, we set cells.z to be negative values.
    ax = get(gcf, 'currentaxes');
    drawContours3(ax, Gt_tmp, -Gt_tmp.cells.z, 20, 'color', 'k');

    % Adjust colorbar
    cmax = max(CO2sat);
    set(findobj(gcf,'type','axes'),'clim',[0, cmax]);

    % Add injection location which is covered up by CO2sat:
    % plot actual location
    plot3(opt.wellXcoord, opt.wellYcoord, -200, 'o', ...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','r',...
        'MarkerSize',10)
    % simulated location
    plot3(var.wellCoord_x, var.wellCoord_y, -200, 'x', ...
        'LineWidth',3,  ...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','k',...
        'MarkerSize',10)

end
    
% -------------------------------------------------------------------------

function plotHammerfestLayers3D(Gt_st, Gt_nd, Gt_tu)
% This plotting routine may be adjusted for desired look:

    G_st = Gt_st.parent;
    G_nd = Gt_nd.parent;
    G_tu = Gt_tu.parent;
    
    figure; set(gcf,'Position',[1 1 1000 800])
    
    % top surface grid of Sto
    %plotGrid(Gt_st,'FaceColor','none','EdgeColor','y')
    
    % three formation grids, colored differently
    %plotGrid(G_st,'FaceColor',[1 .9 .9],'EdgeColor','none')
    hp1 = plotGrid(G_tu,'FaceColor','c','EdgeColor','none');
    hp2 = plotGrid(G_nd,'FaceColor','m','EdgeColor','none');
    hp3 = plotGrid(G_st,'FaceColor',[1 .9 .9],'EdgeColor','none');

    light('Position',[-1 -1 1],'Style','infinite'); lighting phong
    view([-100 30]); grid; axis equal tight off
    hl = legend([hp3,hp2,hp1],{'St\o','Nordmela',['Tub',char(229),'en']},'Location','West');
    set(hl,'Fontsize',20)
    set(gca,'DataAspect',[1 1 0.05]);
    title('Hammerfest Basin Aquifer')
    set(gca,'Fontsize',20)
    
end

% -------------------------------------------------------------------------

function plotHammerfestRockProps(Gt_st, rock2D_st, Gt_nd, rock2D_nd, Gt_tu, rock2D_tu)
    % The following figure can be compared with the plots shown in Compiled
    % CO2 Atlas, chp 6, pg 129.

    
    % Get colorbar limits of the data, for consistency in plotting
    [ntg_min, ntg_max]      = getDataLimits(rock2D_st.ntg, rock2D_nd.ntg, rock2D_tu.ntg);
    [poro_min, poro_max]    = getDataLimits(rock2D_st.poro, rock2D_nd.poro, rock2D_tu.poro);
    [perm_min, perm_max]    = getDataLimits(rock2D_st.perm, rock2D_nd.perm, rock2D_tu.perm);
    [depth_min, depth_max]  = getDataLimits(Gt_st.cells.z, Gt_nd.cells.z, Gt_tu.cells.z);

    
    figure; set(gcf,'Position',[1 1 1480 1011])

    mySubPlot = @(myGt, myData) plotCellData(myGt, myData, 'EdgeColor','none');

    % Sto formation
    subplot(3,4,1)
    mySubPlot(Gt_st, Gt_st.cells.z); colorbar; caxis([depth_min depth_max])
    ylabel('St\o','Fontweight','bold')
    title({'Top surface depth (meter)';' '})
    subplot(3,4,2)
    mySubPlot(Gt_st, rock2D_st.ntg); colorbar; caxis([ntg_min ntg_max]); axis off;
    title({'net-to-gross';' '})
    subplot(3,4,3)
    mySubPlot(Gt_st, rock2D_st.poro); hcbp = colorbar; hcbp.Tag = 'por'; caxis([poro_min poro_max]); axis off;
    title({'porosity';' '})
    subplot(3,4,4)
    mySubPlot(Gt_st, rock2D_st.perm./(milli*darcy)); colorbar; caxis([perm_min perm_max]./(milli*darcy)); axis off;
    title({'permeability (mD)';' '})

    % Nordmela
    subplot(3,4,5)
    mySubPlot(Gt_nd, Gt_nd.cells.z); colorbar; caxis([depth_min depth_max])
    ylabel('Nordmela','Fontweight','bold')
    subplot(3,4,6)
    mySubPlot(Gt_nd, rock2D_nd.ntg); colorbar; caxis([ntg_min ntg_max]); axis off;
    subplot(3,4,7)
    mySubPlot(Gt_nd, rock2D_nd.poro); hcbp = colorbar; hcbp.Tag = 'por'; caxis([poro_min poro_max]); axis off;
    subplot(3,4,8)
    mySubPlot(Gt_nd, rock2D_nd.perm./(milli*darcy)); colorbar; caxis([perm_min perm_max]./(milli*darcy)); axis off;

    % Tubaen
    subplot(3,4,9)
    mySubPlot(Gt_tu, Gt_tu.cells.z); colorbar; caxis([depth_min depth_max])
    ylabel(['Tub',char(229),'en'],'Fontweight','bold')
    subplot(3,4,10)
    mySubPlot(Gt_tu, rock2D_tu.ntg); colorbar; caxis([ntg_min ntg_max]); axis off;
    subplot(3,4,11)
    mySubPlot(Gt_tu, rock2D_tu.poro); hcbp = colorbar; hcbp.Tag = 'por'; caxis([poro_min poro_max]); axis off;
    subplot(3,4,12)
    mySubPlot(Gt_tu, rock2D_tu.perm./(milli*darcy)); colorbar; caxis([perm_min perm_max]./(milli*darcy)); axis off;

    
    % Adjust axis, fonts
    hfig = gcf;
    set(findobj(hfig.Children,'Type','axes'),'Fontsize',16,'box','on')
    axis(findobj(hfig.Children,'Type','axes'),'equal','tight')
    set(findobj(hfig.Children,'Type','axes'),'XTick',[900000 940000 980000],'YTick',[7920000 7960000 8000000 8040000])

    % Adjust ticks in colorbars
    hcbs = findobj(hfig.Children,'Type','colorbar');
    for i = 1:numel(hcbs)
        currTicks = hcbs(i).Ticks;
        set(hcbs(i),'Ticks',[hcbs(i).Limits(1) currTicks hcbs(i).Limits(2)])
        cbt = get(hcbs(i), 'Ticks');
        if strcmpi(hcbs(i).Tag, 'por')
            cbtl = arrayfun(@(t) sprintf('%.2f', t), cbt, 'UniformOutput', false);
        else
            cbtl = arrayfun(@(t) sprintf('%.1f', t), cbt, 'UniformOutput', false);
        end
        set(hcbs(i), 'TickLabels', cbtl)
    end

end

function [dmin, dmax] = getDataLimits(fm1data, fm2data, fm3data)

    dmax = max([    max(fm1data(isfinite(fm1data))), ...
                    max(fm2data(isfinite(fm2data))), ...
                    max(fm3data(isfinite(fm3data)))         ]);

    dmin = min([    min(fm1data(isfinite(fm1data))), ...
                    min(fm2data(isfinite(fm2data))), ...
                    min(fm3data(isfinite(fm3data)))         ]);

end