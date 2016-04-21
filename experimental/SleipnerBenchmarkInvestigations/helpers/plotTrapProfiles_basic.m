function [ hfig ] = plotTrapProfiles_basic( Gt, rock2D, rhoCref, ...
        water_density, seafloor_temp, seafloor_depth, temp_gradient, ...
        press_deviation, sr, sw, dis_max, varargin )

    opt.ZoomIntoIEAGHGregion = true;
    opt.plotAxis = false;
    opt.numContours = 30;
    opt = merge_options(opt, varargin{:});
    
    ta = trapAnalysis(Gt, false);

    % Condensed plot:
    figure; set(gcf,'Position',[1 1 1100 700])
    hfig = gcf;
    

    trapcells = ta.traps~=0;
    cellsTrapVol = zeros(Gt.cells.num,1);
    ta_volumes = volumesOfTraps(Gt, ta);
    cellsTrapVol(trapcells) = ta_volumes(ta.traps(trapcells));


    % GET TRAPPING BREAKDOWN: structural, residual, dissoluion
    % first, compute theoretical capacity (upper bound):
    [ capacityOutput ] = getTrappingCapacities(Gt, rock2D, ta, ...
        rhoCref, water_density, seafloor_temp, seafloor_depth, ...
        temp_gradient, press_deviation, sr, sw, dis_max);

    % Distributed CO2 mass under structural traps: 
    cellsTrapCO2Mass = zeros(Gt.cells.num,1);
    cellsTrapCO2Mass(trapcells) = capacityOutput.strap_mass_co2(trapcells);

    % Cumulative CO2 mass under structural traps:
    trapcaps = accumarray(ta.traps(trapcells), capacityOutput.strap_mass_co2(trapcells));
    trapcap_tot = zeros(Gt.cells.num,1); %ones(size(ta.traps)) * NaN;
    trapcap_tot(trapcells) = trapcaps(ta.traps(trapcells));


    
    %
    subplot(1,2,2); hsub3 = gca; hfsub3 = gcf;
    hold on
    %plotGrid(G, 'EdgeAlpha', 0.1, 'FaceColor', 'none')
    %plotFaces(Gt, bf, 'EdgeColor','k', 'LineWidth',3);
    plotCellData(Gt, trapcap_tot/1e9, trapcap_tot~=0, 'EdgeColor','none')

    % adjust plot
    set(gca,'DataAspect',[1 1 1/100])
    grid; axis tight;
    if opt.ZoomIntoIEAGHGregion
       xlim([436914 440114]); % previously determined
       ylim([6469150 6475050]);
    end
    box
    
    % add title
    title('Capacity', 'FontSize',20)
    % add colorbar
    [ ~ ] = setColorbarHandle( gcf, 'LabelName', 'Mt CO_2', 'fontSize', 22 );
    % adjust axis fontsize
    if opt.plotAxis
       set(gca,'FontSize',15)
    else
       set(gca,'XTickLabel','','YTickLabel','');
    end


    % PLOT SPILL PATHS AND TOPOLOGY
    subplot(1,2,1); hsub5 = gca; hfsub5 = gcf;
    hold on
%     if opt.ZoomIntoIEAGHGregion && strcmpi(name,'Original')
%         % number of contours to draw for GHGT grid
%         % (IEAGHG grid is drawn with 30 contours, however if we are zooming
%         % into IEAGHG region, we first draw mapPlot with more contours,
%         % since once zoomed in, we will see less contours.)
%         numContours = 50;
%     else
%         % number of contours to draw for IEAGHG grid
%         numContours = 30;
%     end
    mapPlot(gcf, Gt, 'traps', ta.traps, 'rivers', ta.cell_lines, ...
        'maplines',opt.numContours);
    %plotFaces(Gt, bf, 'EdgeColor','k', 'LineWidth',3);
    
    % adjust plot
    grid; axis equal tight;
    if opt.ZoomIntoIEAGHGregion
       xlim([436914 440114]); % previously determined
       ylim([6469150 6475050]);
    end
    box

    % add title
    title({'Topography';'and Traps'}, 'FontSize',20)
    % adjust axis fontsize
    if opt.plotAxis
       set(gca,'FontSize',15)
    else
       set(gca,'XTickLabel','','YTickLabel','');
    end


end

