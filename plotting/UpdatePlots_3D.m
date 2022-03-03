function UpdatePlots_3D(model)
    static      = model.static;    
    fig         = static.fig;
    figTitle    = fig.title;
    figStyle    = fig.style;
    figTag      = fig.tag;
    figSubtitle = fig.subtitle;
    f = findobj('type','figure','name',figTitle);
    if(isempty(f))
        f = figure('name',        figTitle, ...
                   'tag',         figTag,   ...
                   'numberTitle', 'off',     ...
                   'windowStyle', figStyle );
        Subplots(f,figSubtitle);
    else
        figure(f);
    end
    if(strcmp(figStyle,'normal')), f.WindowState = 'maximized'; end    
    
    a = findobj(gcf,'type','axes');
    for i = 1:length(a)
        if(strcmp(a(i).Title.String,strcat(''))), set(a(i),'Visible','off'); end
        if(strcmp(a(i).Title.String,strcat(figSubtitle{1}))), ax(1) = a(i); end
        if(strcmp(a(i).Title.String,strcat(figSubtitle{2}))), ax(2) = a(i); end
        if(strcmp(a(i).Title.String,strcat(figSubtitle{3}))), ax(3) = a(i); end
        if(strcmp(a(i).Title.String,strcat(figSubtitle{4}))), ax(4) = a(i); end
        if(strcmp(a(i).Title.String,strcat(figSubtitle{5}))), ax(5) = a(i); end
        if(strcmp(a(i).Title.String,strcat(figSubtitle{6}))), ax(6) = a(i); end
    end  
    %----------------------------------------
    G       = model.grid.G;
    dynamic = model.dynamic;
    params  = dynamic.params;
    qw_inj  = params.qinj(:,1); if(qw_inj(2) ~= 0), qw_inj(1) = qw_inj(2); end
    qo_inj  = params.qinj(:,2); if(qo_inj(2) ~= 0), qo_inj(1) = qo_inj(2); end
    fromIdx = params.fromIdx;
    toIdx   = params.toIdx;
    time    = params.cumScheduleSteps;
    PVI     = params.PVI(1:toIdx);
    pDiff   = params.pDiff(1:toIdx);
    qw      = qw_inj(1:toIdx);
    qo      = qo_inj(1:toIdx);
    Qprod   = params.Qprod;
    Qp_net  = params.Qp_net;
    SwAvg   = params.SwAvg;
    totalPeriod = params.periodEnd(end);
    currentTime = time(toIdx);
    simulation = model.simulation;
    if toIdx == 1
        state = dynamic.states{1,1};
    else
        state = dynamic.states{toIdx-1,1};
    end
    coreLength = model.experiment.geometry.length.value;
    if (isfield(simulation,'bCells'))      
        inlet_mask = G.inlet_mask;
        outlet_mask = G.outlet_mask;
        inlet_frontcells_mask = G.inlet_frontcells_mask;
        outlet_backcells_mask = G.outlet_backcells_mask;
        waterSat = state.s(:,1);
        waterSat(inlet_mask) = waterSat(inlet_frontcells_mask);
        waterSat(outlet_mask) = waterSat(outlet_backcells_mask);
        sw_averaged_over_each_slice = [];
        for i = 1: G.cartDims(1,1)
            mask = G.cells.centroids(:,1) == G.cells.centroids(i,1);
            sw_averaged_over_each_slice = [sw_averaged_over_each_slice;...
                mean(waterSat(mask))];
        end
        x = linspace(0,coreLength,G.cartDims(1,1));
    else
        x = linspace(0,coreLength,G.cartDims(1,1));
        waterSat = state.s(:,1);
        sw_averaged_over_each_slice = [];
        for i = 1: G.cartDims(1,1)
            mask = G.cells.centroids(:,1) == G.cells.centroids(i,1);
            sw_averaged_over_each_slice = [sw_averaged_over_each_slice;...
                mean(waterSat(mask))];
        end
    end
    process     = model.experiment.process;
    processName = lower(process.name);                    
    %----------------------------------------
    displayUnits  = DisplayUnits(model);    
    displaySat    = displayUnits.displaySat;
    displayTime   = displayUnits.displayTime;
    displayLength = displayUnits.displayLength;
    displayPress  = displayUnits.displayPress;
    displayRate   = displayUnits.displayRate;
    displayVolume = displayUnits.displayVolume;
    %----------------------------------------
    satFac = Convert(displaySat);
    waterSat = waterSat ./ satFac;
    timeFac = Convert(displayTime);
    time = time ./ timeFac;
    totalPeriod = totalPeriod ./ timeFac;
    lenFac = Convert(displayLength);
    x = x ./ lenFac;
    rateFrac = Convert(displayRate);
    qo = qo ./ rateFrac;
    qw = qw ./ rateFrac;
    pressFac = Convert(displayPress);
    pDiff = pDiff ./ pressFac;
    volFac = Convert(displayVolume);
    Qprod = Qprod ./ volFac;
    %----------------------------------------                
    updateAxes = model.dynamic.params.updateAxes;
    for k = 1 : length(updateAxes)
        j = updateAxes{k};
        subfigure = figSubtitle{j};
        delete(get(ax(j),'children'))
        set(f,'CurrentAxes',ax(j))
        switch subfigure
            %----------------------------------------
            % visualize saturation map
            %----------------------------------------
            case 'Flooding Experiment' 
                xlabel('x'); ylabel('y'); zlabel('z');
                delete(findall(gcf,'type','annotation'))               
                floodPlot = plotCellData(G,waterSat); view(-39,21); axis equal;    
                h = colorbar;
                t = get(h,'ticks');
                tl = arrayfun(@(x) sprintf('%.2f',x),t,'un',0);
                set(h,'ticklabels',tl)
                if (max(waterSat) > min(waterSat))
                    caxis([min(waterSat) max(waterSat)]);
                end
                if(isfield(model.plot,'colormap'))
                    if(strcmp(model.plot.colormap,'jet')), colormap(jet); end                
                    if(strcmp(model.plot.colormap,'parula')), colormap(parula(64)); end
                    if(strcmp(model.plot.colormap,'hsv')), colormap(hsv); end
%                     if(strcmp(processName,'drainage')), colormap(flipud(colormap)); end
%                     if(strcmp(processName,'bl')), colormap(flipud(colormap)); end
                end
                caxis([0 1])
                dim = [.05 0.5 0.3 0.4];
                caxis([0 1]);
                duration = seconds(currentTime);
                duration.Format = 'hh:mm:ss';
                str1 = strcat('Time', {' = '}, cellstr(duration));
                format short g; a = sprintf('%.2f', PVI(end)); b = str2double(a);
                str2 = strcat('PVI =', {' '}, num2str(b));
                str = [str1;str2];
                annotation('textbox',dim,'string',str,'edgecolor','none','fontsize',9); 
                getframe(gcf);
            %----------------------------------------
            % visualize saturation front
            %----------------------------------------
            case 'Saturation Front'                
%                 plotToolbar(model.grid.G, dynamic.states, 'field', 's:1', 'plot1d', true, ...
%                            'lockCaxis', true);
                marker  = 'b-';
                markerSize = 4;
                lineWidth = 1.5;
                satFrontPlot = plot(x,sw_averaged_over_each_slice,marker, 'MarkerSize',markerSize,...
                                    'MarkerFaceColor','r',...
                                    'LineWidth',lineWidth);
                xlim([0 max(x)])
                ylim([min(0) max(1)])                             
                grid on
                title(subfigure);
                xstr = strcat('$$ \it {\bf {Distance}}',{' '},'[',displayLength,'] $$');
                xlabel(xstr,'interpreter','latex')
                ystr = strcat('$$ \it {\bf {S_{w}}} $$');
                ylabel(ystr,'Interpreter','latex');
                yTicks = (0:0.2:1);
                set(gca,'YTick',yTicks)
                set(gca,'XMinorTick','on')
        end
    end
    
    selectedTimeAxes = model.dynamic.params.selectedTimeAxes;
    for k = 1 : length(selectedTimeAxes)
        j = selectedTimeAxes{k};
        subfigure = figSubtitle{j};
        set(f,'CurrentAxes',ax(j))
        lineWidth = 1.5;
        lineStyle = '-';
        lineColor = [0.5 0.5 0.5];
        labelVerticalAlignment = 'middle';
        labelHorizontalAlignment = 'center';
        duration = seconds(currentTime);
        duration.Format = 'hh:mm:ss';
        label = strcat(cellstr(duration));
        fontSize = 7;
        xValue = time(toIdx);
        switch subfigure
            %----------------------------------------
            % current time injection rates
            %----------------------------------------
            case 'Injection Rates'
                injLine = findobj(ax(j),'type','constantline');
                if ~isempty(injLine), delete(injLine); end
                injLine = xline(xValue,lineStyle,label, ...
                    'linewidth',lineWidth, ...
                    'fontsize',fontSize, ...
                    'labelverticalalignment',labelVerticalAlignment, ...
                    'labelhorizontalalignment',labelHorizontalAlignment, ...
                    'color',lineColor, ...
                    'tag','injLine');
            %----------------------------------------
            % current time pressure differential
            %----------------------------------------
            case 'Pressure Differential'
                pressLine = findobj(ax(j),'type','constantline');
                if ~isempty(pressLine), delete(pressLine); end
                pressLine = xline(xValue,lineStyle,label, ...
                    'linewidth',lineWidth, ...
                    'fontsize',fontSize, ...
                    'labelverticalalignment',labelVerticalAlignment, ...
                    'labelhorizontalalignment',labelHorizontalAlignment, ...
                    'color',lineColor, ...
                    'tag','pressLine');
            %----------------------------------------
            % current time average water saturation
            %----------------------------------------
            case 'Average Water Saturation'
                swAvgLine = findobj(ax(j),'type','constantline');
                if ~isempty(swAvgLine), delete(swAvgLine); end
                swAvgLine = xline(xValue,lineStyle,label, ...
                    'linewidth',lineWidth, ...
                    'fontsize',fontSize, ...
                    'labelverticalalignment',labelVerticalAlignment, ...
                    'labelhorizontalalignment',labelHorizontalAlignment, ...
                    'color',lineColor, ...
                    'tag','swAvgLine');               
            %----------------------------------------
            % current time Cummulative Production, Np and NpEff
            %----------------------------------------
            case 'Cummulative Production'
                cumLine = findobj(ax(j),'type','constantline');
                if ~isempty(cumLine), delete(cumLine); end
                cumLine = xline(xValue,lineStyle,label, ...
                    'linewidth',lineWidth, ...
                    'fontsize',fontSize, ...
                    'labelverticalalignment',labelVerticalAlignment, ...
                    'labelhorizontalalignment',labelHorizontalAlignment, ...
                    'color',lineColor, ...
                    'tag','cumLine');
                legend('off')
        end
    end
end