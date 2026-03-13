function filtered_data = PlotObservation_pre_hm(model, filter_no)
%
% DESCRIPTION: plots the experimental measurements before the history
%              matching starts, and applies also a median filter to it,
%              which can be replaced with the original data to improve the
%              efficiency of the history matching
%
% SYNOPSIS:
%   filtered_data = PlotObservation_pre_hm(model, filter_no)
%
% PARAMETERS:
%   model - struct containing following fields:
%   - experiment
%
% RETURNS:
%    - plot of the experimental measurements with the median filter on the
%    top
%   - filtered_data: the results of the median filter on the data
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
    static = SetupStaticParams(model); % set static struct   
    fig = static.fig;
    figTitle = fig.title;    
    figSubtitle = fig.subtitle;
    
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
    f = figure('Name', 'Exploratory data analysis','NumberTitle', 'off');
    colormap(f, gray)
    tiledlayout('flow')
    %----------------------------------------

    marker  = 'ks';
    markerSize = 4;
    lineWidth = 0.5;
    
    if(isfield(observation,'pressure'))
        if(observation.pressure.include)
            nexttile;
            time_pdif = observation.pressure.table{:,1};
            time_pdif = time_pdif ./ timeFac;            
            pDiff = observation.pressure.table{:,2};
            pDiff = pDiff ./ pressFac;
            %----------------------------------------
            % observed pressure differential data
            %----------------------------------------
            plot(time_pdif, pDiff,marker,'MarkerSize',markerSize,'LineWidth',lineWidth);                               
            grid off; title('Pressure Differential');
            xlabel('Time (' + string(displayTime) + ')')
            ylabel('Pressure (' + string(displayPress) + ')')

            hold on
            filtered_pressure = medfilt1(pDiff, filter_no.pressure);
            % save the data in the same unit as input
            filtered_data.pressure = filtered_pressure * pressFac;
            plot(time_pdif, filtered_pressure, 'LineWidth',2)
            legend('Original data', 'Median Filter', 'Location','best')

        end
    end
    
    if(isfield(observation,'swavg'))
        if(observation.swavg.include)
            nexttile;
            time_swavg = observation.swavg.table{:,1};
            time_swavg = time_swavg ./ timeFac;
            SwAvg = observation.swavg.table{:,2};
            SwAvg = SwAvg ./ satFac;
            %----------------------------------------
            % observed average water saturation data
            %----------------------------------------
            plot(time_swavg,SwAvg,marker,'MarkerSize',markerSize,'LineWidth',lineWidth);
            grid off; title('Average Water Saturation')
            xlabel('Time (' + string(displayTime) + ')')
            ylabel('Water saturation (' + string(displaySat) + ')')

            hold on
            filtered_sw_avg =medfilt1(SwAvg, filter_no.sw_avg);
            % save the data in the same unit as input
            filtered_data.sw_avg = filtered_sw_avg * satFac;
            plot(time_swavg, filtered_sw_avg, 'LineWidth',2)
            legend('Original data', 'Median Filter', 'Location','best')
        end
    end
    
    if(isfield(observation,'prod'))
        if(observation.prod.include)
            nexttile;
            time_prod = observation.prod.table{:,1};
            time_prod = time_prod ./ timeFac;
            Qprod = observation.prod.table{:,2};
            Qprod = Qprod ./ volFac;
            %----------------------------------------
            % observed cumulative production data
            %----------------------------------------          
            plot(time_prod, Qprod,marker,'MarkerSize',markerSize,'LineWidth',lineWidth);
            grid off; title('Production')
            xlabel('Time (' + string(displayTime) + ')')
            ylabel('Production (' + string(displayVolume) + ')')

            hold on
            filtered_prod = medfilt1(Qprod, filter_no.prod);
            % save the data in the same unit as input
            filtered_data.prod = filtered_prod * volFac;
            plot(time_prod, filtered_prod, 'LineWidth',2)
            legend('Original data', 'Median Filter', 'Location','best')

        end
    end
end