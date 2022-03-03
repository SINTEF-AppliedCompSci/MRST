function Visualize_3D(model)
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
    state = dynamic.states{toIdx-1,1};
    simulation = model.simulation;
    coreLength = model.experiment.geometry.length.value;
    observation = model.experiment.observation;
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
    Qp_net = Qp_net ./ volFac;    
    %----------------------------------------
%     if isfield(model.experiment,'observation')
%         if(~isempty(model.experiment.observation))
%             PlotObservation(model);
%         end     
%     end
    %----------------------------------------
    for j = 1 : length(figSubtitle)
        subfigure = figSubtitle{j};
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
%                     if(strcmp(processName,'imbibition')), colormap(flipud(colormap)); end
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
            % visualize injection rates
            %----------------------------------------
            case 'Injection Rates'
                markerSize = 4;
                lineWidth = 1.5;               
                minY = min([qo;qw]);
                maxY = max([qo;qw]); 
                axRate = gca;
                for i = fromIdx : toIdx - 1
                    marker  = 'r-';
                    qoPlot = plot(time(i:i+1),qo(i:i+1),marker,...
                                  'MarkerSize',markerSize,...
                                  'LineWidth',lineWidth);
                    xlim([0 totalPeriod])
                    if(axRate.YLim(2) < 1.1 * maxY), ylim([-.1 1.1 * maxY]); end                   
                    hold on % qw
                    marker  = 'b-';
                    qwPlot = plot(time(i:i+1),qw(i:i+1),marker, ...
                                  'MarkerSize',markerSize,...
                                  'LineWidth',lineWidth);
                end
                grid on
                title(subfigure)
                xstr = strcat('$$ \it {\bf {Time}}',{' '},'[',displayTime,'] $$');
                xlabel(xstr,'interpreter','latex')
                ystr = strcat('$$ \it {\bf {q}}',{' '},'[',displayRate,'] $$');
                ylabel(ystr,'Interpreter','latex','Color','k');                
                set(gca,'XMinorTick','on')
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
            %----------------------------------------
            % visualize pressure differential
            %----------------------------------------
            case 'Pressure Differential'                 
                marker  = 'bo-';
                markerSize = 4;
                lineWidth = 1.5;                            
                minY = min(pDiff);
                maxY = max(pDiff);
                axPress = gca;
                for i = fromIdx : toIdx - 1
                    pressPlot = plot(time(i:i+1),pDiff(i:i+1),marker,...
                                     'MarkerSize',markerSize,...
                                     'LineWidth',lineWidth);
                    hold on
                end                
                xlim([0 totalPeriod])
                if(axPress.YLim(2) < 1.1 * maxY), ylim([0.9 * minY 1.1 * maxY]); end
                grid on
                title(subfigure)
                xstr = strcat('$$ \it {\bf {Time}}',{' '},'[',displayTime,'] $$');
                xlabel(xstr,'interpreter','latex')
                ystr = strcat('$$ \it {\bf {\Delta P}}',{' '},'[',displayPress,'] $$');
                ylabel(ystr,'Interpreter','latex');
                set(gca,'XMinorTick','on')
            %----------------------------------------
            % visualize average water saturation
            %----------------------------------------
            case 'Average Water Saturation'
                marker  = 'bo-';
                markerSize = 4;
                lineWidth = 1.5;
                for i = fromIdx : toIdx - 1
                    satAvgPlot = plot(time(i:i+1),SwAvg(i:i+1),marker,...
                                        'MarkerSize',markerSize,...
                                        'LineWidth',lineWidth);
                    hold on
                end
                xlim([0 totalPeriod])
                ylim([0 1])
                grid on
                title(subfigure)
                xstr = strcat('$$ \it {\bf {Time}}',{' '},'[',displayTime,'] $$');
                xlabel(xstr,'interpreter','latex')
                ystr = strcat('$$ \it {\bf {\overline {S_{w}}}} $$');
                ylabel(ystr,'Interpreter','latex');
                yTicks = (0:0.2:1);
                set(gca,'YTick',yTicks)
                set(gca,'XMinorTick','on')
            %----------------------------------------
            % visualize Cummulative Production, Np and NpEff
            %----------------------------------------            
            case 'Cummulative Production'           
                if(processName == "drainage")
                    % water production
                    Q = Qprod(:,1); Qnet = Qp_net(:,1);
                    markerProd  = 'b-o'; markerNet = 'k-o'; 
                    axColorProd = 'b'; axColorNet = 'k';
                end
                if(processName == "imbibition")
                    % oil production
                    Q = Qprod(:,2); Qnet = Qp_net(:,2);
                    markerProd  = 'b-o'; markerNet = 'k-o'; 
                    axColorProd = 'b'; axColorNet = 'k';
                end
                markerSize = 4;
                lineWidth = 1.5;
                hLegend = findobj(f,'tag','legend');
                if(~isempty(hLegend)), legend('visible','off'); end
                for i = fromIdx : toIdx - 1                   
                    yyaxis right
                    minY = min(Q); maxY = max(Q); 
                    QPlot = plot(time(i:i+1),Q(i:i+1),markerProd, ...
                                  'MarkerSize',markerSize,...
                                  'LineWidth',lineWidth);                    
                    ax2 = gca; % current axes
                    ax2.YColor = axColorProd;
                    xlim([0 totalPeriod])
                    if(ax2.YLim(2) < 1.1 * maxY), ylim([-.1 1.1 * maxY]); end
                    if (maxY > 0), ylim([0.9 * minY 1.1 * maxY]); end
                    
                    yyaxis left
                    minY = min(Qnet); maxY = max(Qnet);
                    QnetPlot = plot(time(i:i+1),Qnet(i:i+1),markerNet,...
                                  'MarkerSize',markerSize,...
                                  'LineWidth',lineWidth);
                    ax1 = gca; % current axes
                    ax1.YColor = axColorNet;                    
                    if(ax1.YLim(2) < 1.1 * maxY), ylim([-.1 1.1 * maxY]); end
                    if (maxY > 0), ylim([0.9*minY 1.1*maxY]); end                    
                    hold on        
                end
                if(isempty(hLegend))
                    legend([QnetPlot QPlot],{'Production-Injection','Production'},'location','southeast')
                else
                    if(isfield(observation,'prod'))
                        if(observation.prod.include)
                            hLegend.String(4:end) = [];                        
                        else
                            hLegend.String(3:end) = [];
                        end
                    else
                        hLegend.String(3:end) = [];
                    end
                end
                legend('show')
                grid on
                title(subfigure)
                xstr = strcat('$$ \it {\bf {Time}}',{' '},'[',displayTime,'] $$');
                xlabel(xstr,'interpreter','latex')
                ystr = strcat('$$ \it {\bf {Q}}',{' '},'[',displayVolume,'] $$');
                ylabel(ystr,'Interpreter','latex','Color','k');                
                set(gca,'XMinorTick','on')
            otherwise
                disp('axis not found')
        end        
    end    
%     if(strcmp(figStyle,'normal')), f.WindowState = 'minimized'; end

    if isfield(model.experiment,'observation')
        if(~isempty(model.experiment.observation))
            PlotObservation(model);
        end     
    end
    ShowSlider_3D(model);
end