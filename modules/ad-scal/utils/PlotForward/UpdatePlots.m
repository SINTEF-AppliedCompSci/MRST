function UpdatePlots(model)
%
% DESCRIPTION: updates the plots created with the Visualize module when the
% slider is moved a time step is chosen from the list in the figure
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
%   updates the simulation prediction plots based on user interaction with
%   slider or time step lists in the simulation prediction plots
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
    qw      = qw_inj(1:toIdx);
    qo      = qo_inj(1:toIdx);
    Qprod   = params.Qprod;
    Qp_net  = params.Qp_net;
    SwAvg   = params.SwAvg;
    totalPeriod = params.periodEnd(end);
    currentTime = time(toIdx);
    simulation = model.simulation;
    if toIdx == 1
        state = model.state0;
    else
        state = dynamic.states{toIdx-1,1};
    end
    coreLength = model.experiment.geometry.length.value;
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
                floodPlot = plotCellData(G,waterSat); axis equal;    
                h = colorbar;
                t = get(h,'ticks');
                tl = arrayfun(@(x) sprintf('%.2f',x),t,'un',0);
                set(h,'ticklabels',tl)
                dim = [.05 0.5 0.3 0.4];
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
                ylim([0 1])                             
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
                        hold on
                        if any(observed_sat_profile_t ./ timeFac == time(toIdx))
                            satFrontPlot_observation = plot(observed_sat_profile_location,...
                                observed_sat_profile(observed_sat_profile_t ./ timeFac == time(toIdx),:)...
                                ,'ro', 'MarkerSize',2);
                        end
                    end
                end
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
            % current time Production, Np and NpEff
            %----------------------------------------
            case 'Production'
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