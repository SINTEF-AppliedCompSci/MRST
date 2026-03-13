function PlotObservation(model)
%
% DESCRIPTION: adds the experimental measurements to the plots created by
% the Visualize module
%
% SYNOPSIS:
%   PlotObservation(model)
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
%   experimental measurements added to the simulation prediction plots
%   created with the Visualize modules
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
    static = model.static;    
    fig = static.fig;
    figTitle = fig.title;    
    figSubtitle = fig.subtitle;
    f = findobj('type','figure','name',figTitle);
    
    % display units
    displayUnits  = DisplayUnits(model);    
    displaySat    = displayUnits.displaySat;
    displayTime   = displayUnits.displayTime;
    displayLength = displayUnits.displayLength;
    displayPress  = displayUnits.displayPress;
    displayRate   = displayUnits.displayRate;
    displayVolume = displayUnits.displayVolume;

    % unit converter
    satFac   = Convert(displaySat);
    timeFac  = Convert(displayTime);
    pressFac = Convert(displayPress);
    lenFac   = Convert(displayLength);
    rateFrac = Convert(displayRate);
    volFac   = Convert(displayVolume);

    %----------------------------------------
    observation = model.experiment.observation;
    figure(f);
    a = findobj(f,'type','axes');
    for i = 1:length(a)
        if(strcmp(a(i).Title.String,strcat(''))), set(a(i),'Visible','off'); end
        if(strcmp(a(i).Title.String,strcat('Pressure Differential'))), axPress = a(i); end % 'Pressure Differential'
        if(strcmp(a(i).Title.String,strcat('Average Water Saturation'))), axSwAvg = a(i); end % 'Average Water Saturation'
        if(strcmp(a(i).Title.String,strcat('Production'))), axQprod = a(i); end % 'Production'
    end 
    %----------------------------------------
    
    params   = model.dynamic.params;
    fromTime = params.periodStart(end) / timeFac;
    toTime   = params.periodEnd(end) / timeFac;
    
    process     = model.experiment.process; 
    
    if(isfield(observation,'pressure'))
        if(observation.pressure.include)            
            time_pdif = observation.pressure.table{:,1};
            time_pdif = time_pdif ./ timeFac;            
            pDiff = observation.pressure.table{:,2};
            pDiff = pDiff ./ pressFac;
            if isfield(model.experiment.observation, "pressure_mid")
                time_pdif_mid = observation.pressure_mid.table{:,1};
                time_pdif_mid = time_pdif_mid ./ timeFac; 
                pDiff_mid = observation.pressure_mid.table{:,2};
                pDiff_mid = pDiff_mid ./ pressFac;
            end
            if ishghandle(axPress) % pressure results exist
                %----------------------------------------
                % observed pressure differential data
                %----------------------------------------
                set(f,'CurrentAxes',axPress)    
                marker  = 'ks';
                markerSize = 2;
                lineWidth = 0.5;
                xLim = axPress.XLim;
                yLim = axPress.YLim;
                idx = find(time_pdif >= fromTime & time_pdif < toTime);
                plot(time_pdif(idx),pDiff(idx),marker,'MarkerSize',markerSize,'LineWidth',lineWidth); 
                if isfield(model.experiment.observation, "pressure_mid")
                idx_Pmid = find(time_pdif_mid >= fromTime & time_pdif_mid < toTime);
                plot(time_pdif_mid(idx_Pmid),pDiff_mid(idx_Pmid),marker,'MarkerSize',markerSize,'LineWidth',lineWidth); 
                end
                xlim(xLim);
                ylim(yLim);
                if(yLim(2) < max(pDiff)), ylim([yLim(1) 1.1 * max(pDiff)]); end
            end
        end
    end
    
    if(isfield(observation,'swavg'))
        if(observation.swavg.include)
            time_swavg = observation.swavg.table{:,1};
            time_swavg = time_swavg ./ timeFac;
            SwAvg = observation.swavg.table{:,2};
            SwAvg = SwAvg ./ satFac;
            if ishghandle(axSwAvg) % saturation results exist
                %----------------------------------------
                % observed average water saturation data
                %----------------------------------------
                set(f,'CurrentAxes',axSwAvg)
                marker  = 'ro';
                markerSize = 2;
                lineWidth = 0.5;
                xLim = axSwAvg.XLim;
                yLim = axSwAvg.YLim;
                idx = find(time_swavg >= fromTime & time_swavg < toTime);
                plot(time_swavg(idx),SwAvg(idx),marker,'MarkerSize',markerSize,'LineWidth',lineWidth);
                xlim(xLim);
                ylim(yLim);
                grid on;
            end
        end
    end
    
    if(isfield(observation,'prod'))
        if(observation.prod.include)
            time_prod = observation.prod.table{:,1};
            time_prod = time_prod ./ timeFac;
            Qprod = observation.prod.table{:,2};
            Qprod = Qprod ./ volFac;
            if ishghandle(axQprod) % production results exist
                %----------------------------------------
                % observed cumulative production data
                %----------------------------------------
                set(f,'CurrentAxes',axQprod)            
                marker  = 'ks';
                markerSize = 4;
                lineWidth = 0.5;
                xLim = axQprod.XLim;
                yLim = axQprod.YLim;
                idx = find(time_prod >= fromTime & time_prod < toTime);
                plot(time_prod(idx),Qprod(idx),marker,'MarkerSize',markerSize,'LineWidth',lineWidth);
                
                if not(strcmpi(process.type, 'cent'))
                    
                    xlim(xLim);
                    if axQprod.YLim(2) > max(Qprod)
                        ylim(yLim);
                    else
                        ylim([0 max(Qprod)]);
                    end
                    
                    hLegend = findobj(f,'tag','legend');
                    if(length(hLegend.String) <= 3)
                        hLegend.String{3} = 'Obs';                    
                    else
                        hLegend.String(4:end) = [];
                    end
                    if(length(hLegend.String) <= 2)
                        hLegend.String{2} = 'Experimental Measurement';                    
                    else
                        hLegend.String(3:end) = [];
                    end
                    
                else
                    
                    hLegend = findobj(f,'tag','legend');
                    if(length(hLegend.String) <= 2)
                        hLegend.String{2} = 'Experimental Measurement';                    
                    else
                        hLegend.String(3:end) = [];
                    end
                end
            end
        end
    end
end