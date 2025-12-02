function Visualize(model)
%
% DESCRIPTION: Visualization of the simulation results after each schedule
%              row
%
% SYNOPSIS:
%   Visualize(model)
%
% PARAMETERS:
%   model - struct containing following fields:
%   - state: states of the simulation in MRST format
%   - grid
%   - dynamic: information saved during the simulation
%   - experiment: saturation functions used for forward modeling
%   - simulation: time stepping and grid cells information
%
% RETURNS:
%   plot of the simulation results, including saturation, pressure,
%   production, and rates
%
% ----------------------------------
% (c) 2020-2022
% Siroos Azizmohammadi
% Omidreza Amrollahinasab
% Montanuniversität Leoben, Austria
% Chair of Reservoir Engineering
% https://dpe.ac.at/
% ----------------------------------
%
%%
static      = model.static;    
fig         = static.fig;
figTitle    = fig.title;
if isfield(fig, 'style')
    figStyle    = fig.style;
else
    figStyle = 'docked';
end
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
if isfield(model.experiment.observation, "pressure_mid")
    pDiff_mid = params.pDiff_mid(1:toIdx);
end
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
    waterSat = state.s(:,1);
    waterSat(1) = state.s(2,1);
    waterSat(end) = state.s(end-1,1);
    waterSatFull = [interp1(G.cells.centroids(2:3,1),state.s(2:3,1),0,'linear','extrap');...
        state.s(2:end-1,1);interp1(G.cells.centroids(end-2:end-1,1) - ...
        2 * G.cells.centroids(1,1),state.s(end-2:end-1,1),coreLength,'linear','extrap')];
    x = [0;G.cells.centroids(2:end-1,1) - 2 * G.cells.centroids(1,1);coreLength];
else
    x = [0;G.cells.centroids(:,1);coreLength];
    waterSat = state.s(:,1);
    waterSatFull = [state.s(1,1);state.s(:,1);state.s(end,1)];
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
if isfield(model.experiment.observation, "pressure_mid")
    pDiff_mid = pDiff_mid ./ pressFac;
end
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
            if (params.Sw_max > params.Sw_min)
                caxis([params.Sw_min params.Sw_max]);
            end
            color_map = model.plot.colormap.inputColormap;
            if(isfield(model.plot,'colormap'))
                if(strcmp(color_map,'jet')), colormap(jet); end                
                if(strcmp(color_map,'parula')), colormap(parula(64)); end
                if(strcmp(color_map,'hsv')), colormap(hsv); end
                if(strcmp(color_map,'winter')), colormap(flip(winter)); end
            end
            dim = [.05 0.5 0.3 0.4];
            duration = seconds(currentTime);
            duration.Format = 'hh:mm:ss';
            str1 = strcat('Time', {' = '}, cellstr(duration));
            format short g; a = sprintf('%.2f', PVI(end)); b = str2double(a);
            str2 = strcat('PVI =', {' '}, num2str(b));
            str = [str1;str2];
%             annotation('textbox',dim,'string',str,'edgecolor','none','fontsize',9);
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
%         plotToolbar(model.grid.G, dynamic.states, 'field', 's:1', 'plot1d', true, ...
%                    'lockCaxis', true);
            marker  = 'k-';
            markerSize = 4;
            lineWidth = 1.5;
            if model.experiment.rock.heterogeneous
                satFrontPlot = stairs(x,waterSatFull,marker, 'MarkerSize',markerSize,...
                                    'MarkerFaceColor','r',...
                                    'LineWidth',lineWidth);
            else
                satFrontPlot = plot(x,waterSatFull,marker, 'MarkerSize',markerSize,...
                    'MarkerFaceColor','r',...
                    'LineWidth',lineWidth);
            end
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
            % visualize saturation front observations
            %----------------------------------------   
            if isfield(model.experiment.observation,'satProfile')
                if isfield(model.experiment.observation.satProfile,'table')
                    observed_sat_profile_t = model.experiment.observation.satProfile.table(2:end,1);
                    observed_sat_profile_location = model.experiment.observation.satProfile.table(1,2:end) ./ lenFac;
                    observed_sat_profile = model.experiment.observation.satProfile.table(2:end,2:end);
                    if any(observed_sat_profile_t ./ timeFac == time(end))
                        hold on
                        satFrontPlot_observation = plot(observed_sat_profile_location,...
                            observed_sat_profile(observed_sat_profile_t ./ timeFac == time(end),:)...
                            ,'ro', 'MarkerSize',2); hold off
                    end
                end
            end
        %----------------------------------------
        % visualize pressure differential
        %----------------------------------------
        case 'Pressure Differential'                 
            marker  = 'r-';
            markerSize = 2;
            lineWidth = 1.5;                            
            minY = min(pDiff);
            maxY = max(pDiff);
            axPress = gca;
            for i = fromIdx : toIdx - 1
                pressPlot = plot(time(i:i+1),pDiff(i:i+1),marker,...
                                 'MarkerSize',markerSize,...
                                 'LineWidth',lineWidth);
                if isfield(model.experiment.observation, "pressure_mid")
                press_mid_Plot = plot(time(i:i+1),pDiff_mid(i:i+1),marker,...
                                 'MarkerSize',markerSize,...
                                 'LineWidth',lineWidth);
                end
                hold on
            end                
            xlim([0 totalPeriod])
            if(axPress.YLim(2) < 1.1 * maxY), ylim([0.9 * minY 1.1 * maxY]); end
            grid off
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
            marker  = 'k-o';
            markerSize = 2;
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
        case 'Production'           
            
            if(processName == "drainage")
                % water production
                Q = Qprod(:,1); Qnet = Qp_net(:,1);
                markerProd  = 'r-'; markerNet = 'k-'; 
                axColorProd = 'r'; axColorNet = 'k';
            end
            if(processName == "imbibition")
                % oil production
                Q = Qprod(:,2); Qnet = Qp_net(:,2);
                markerProd  = 'r-'; markerNet = 'k-'; 
                axColorProd = 'r'; axColorNet = 'k';
            end
            markerSize = 2;
            lineWidth = 1.5;
            hLegend = findobj(f,'tag','legend');
            if(~isempty(hLegend)), legend('visible','off'); end
            if not(strcmpi(process.type, 'cent')) && not(strcmpi(process.type, 'uss'))
                for i = fromIdx : toIdx - 1                   
                    yyaxis right
                    minY = min(Q); maxY = max(Q); 
                    QPlot = plot(time(i:i+1),Q(i:i+1),markerProd, ...
                                  'MarkerSize',markerSize,...
                                  'LineWidth',lineWidth);                    
                    ax2 = gca; % current axes
                    ax2.YColor = axColorProd;
                    xlim([0 totalPeriod])
                    if (maxY > 0), ylim([0.9 * minY 1.1 * maxY]); end
                    ystr = strcat('$$ \it {\bf {Q}}',{' '},'[',displayVolume,'] $$');
                    ylabel(ystr,'Interpreter','latex','Color','r');  

                    yyaxis left
                    minY = min(Qnet); maxY = max(Qnet);
                    QnetPlot = plot(time(i:i+1),Qnet(i:i+1),markerNet,...
                                  'MarkerSize',markerSize,...
                                  'LineWidth',lineWidth);
                    ax1 = gca; % current axes
                    ax1.YColor = axColorNet;                    
                    if (maxY > 0), ylim([0.9*minY 1.1*maxY]); end                    
                    hold on        

                end
                if(isempty(hLegend))
                    legend([QnetPlot QPlot],{'Production-Injection',' Cumulative Production'},'location','southeast')
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
            else
                for i = fromIdx : toIdx - 1                   
                    minY = min(Qnet); maxY = max(Qnet);
                    QnetPlot = plot(time(i:i+1),Qnet(i:i+1),markerNet,...
                                  'MarkerSize',markerSize,...
                                  'LineWidth',lineWidth);
                    ax1 = gca; % current axes
                    ax1.YColor = axColorNet;                    
                    if (maxY > 0), ylim([0.9*minY 1.1*maxY]); end  
                    xlim([0 totalPeriod])
                    hold on        
                end
                if(isempty(hLegend))
                    legend(QnetPlot,'Production-Injection','location','southeast')
                else
                    if(isfield(observation,'prod'))
                        if(observation.prod.include)
                            hLegend.String(3:end) = [];                        
                        else
                            hLegend.String(2:end) = [];
                        end
                    else
                        hLegend.String(2:end) = [];
                    end
                end
            end
            legend('show')
            grid off
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
ShowSlider(model);