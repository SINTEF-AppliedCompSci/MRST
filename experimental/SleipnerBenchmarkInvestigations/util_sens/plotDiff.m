function [ plumes_comb ] = plotDiff( Gt, plumes, sim_report, states, topsurface )

    
    %% Combining plume info together for the same year
    % i.e., if plumes{i}.year is the same for more than one index of i,
    % than hCO2 are added together (assuming they correspond to
    % non-overlapping polygons) and polygon outlines are stored into one
    % field called outlines
    
    plume_years = [1999; 2001; 2002; 2004; 2006; 2008];
    plumes_comb = cell(numel(plume_years),1);
    for i = 1:numel(plume_years)

        % pre-allocate with zeros
        hCO2_combined = zeros(Gt.cells.num,1);

        for j = 1:numel(plumes)
            if plumes{j}.year == plume_years(i);
                hCO2_combined = plumes{j}.h + hCO2_combined;
                outlines{j} = plumes{j}.outline;
            end
        end

        plumes_comb{i}.hCO2 = hCO2_combined;
        plumes_comb{i}.year = plume_years(i);
        plumes_comb{i}.outlines = outlines(~cellfun('isempty',outlines));
        clear outlines
    end

    %% Plotting
    % Here, the first subplot shows hCO2 and the polygon outlines. The
    % second subplot shows the simulated CO2 heights under the caprock. The
    % third subplot shows the difference between the heights in the first
    % and second subplots. All subplots coorespond to the same year.
    
    myplotCellData =@(G,data) plotCellData(G,data,'EdgeColor','none');
    for i=1:numel(plumes_comb)

        % hack
        %pl_yr = plumes_comb{i}.year - 1998;
        %tstep = find(sim_year==pl_yr);
        
        
        % get reservoir time index
        [rti,~] = find(sim_report.ReservoirTime/year == (plume_years(i) - plume_years(1) + 1) );
        
        if(rti<=numel(states))
            figure; set(gcf,'Position',[1 1 1500 500]); %figure(i),clf

            subplot(1,3,1)
            hold on
            myplotCellData(Gt, plumes_comb{i}.hCO2); colorbar; axis equal tight;

            for j=1:numel(plumes_comb{i}.outlines)
                onePolygon = plumes_comb{i}.outlines{j};
                line(onePolygon(:,1), onePolygon(:,2), topsurface(onePolygon)-1, 'LineWidth',2, 'Color','r');
            end

            title({'hCO2 = (topfit(x,y) - topsurface(x,y)) > 0, of polygon(s)';['in year ',num2str(plumes_comb{i}.year)]})
            legend('hCO2 of polygon','polygon outline(s)')

            subplot(1,3,2)
            myplotCellData(Gt, states{rti}.h); colorbar; axis equal tight;
            title({'simulated CO2 height under top surface';['in year ',num2str(plumes_comb{i}.year)]})

            subplot(1,3,3);
            myplotCellData(Gt, states{rti}.h - plumes_comb{i}.hCO2); colorbar; axis equal tight;
            title('simulated CO2 height minus hCO2 of polygon(s)')

        end

    end


end

