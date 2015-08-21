function [ hfig, hax ] = makeSideProfilePlots( Years2plot, inj_year, schedule, G, Gt, sim_report, states, rock2D, fluid, rhoCref, wellCellIndex, trapstruct )


% Reservoir Time2plot is an array of all times in seconds you wish to
% visualize the CO2 plume (saturation). If any times coorespond to the
% observation plume data, that plume outline is plotted as well.

ReservoirTime2plot  = (Years2plot - inj_year(1)+1 ).*(365*24*60*60); % seconds

% total mass injected during the schedule (included both injection and
% migration periods)
[ hfig, ~, timeSinceInj, massNow ] = plotInjectRateVsTime(schedule, inj_year, rhoCref); % massNow is in Mt
timeSinceInj_s  = timeSinceInj * year; % in seconds
close(hfig);


%figure; %set(gcf, 'Position', [1 1 1500 1100])
%hold on

for i = 1:numel(ReservoirTime2plot)
    
    % get reservoir time index
    [rti,~] = find(sim_report.ReservoirTime==ReservoirTime2plot(i));
    
    % Note: Each reservoir time gets a separate figure for plotting
    
    % Prepare for Panel Plots:
    sol.h       = zeros(Gt.cells.num, 1);
    [ii, jj]    = ind2sub(G.cartDims, G.cells.indexMap);
    opts = {'slice',     double([ii(wellCellIndex), jj(wellCellIndex)]),...
            'Saxis',     [0 1-fluid.res_water], ...
            'view',      [-85 70],...
            'plotPlume', true, ...
            'wireH',     true,...
            'wireS',     true,...
            'maxH',      (max(Gt.cells.z) - min(Gt.cells.z))/3, ...
            'plotHist',  false};
    t = 0;
    W = schedule.control(1).W;
    
    % TODO: ensure appropriate figure plotted on
    plotPanelVE(G, Gt, W, sol, t, [0,0,0,0,0,0], opts{:});


    % Then make Panel Plot for a given time:
    W           = schedule.control( schedule.step.control(rti) ).W;
    state{:}    = states{ schedule.step.control(rti) };
    t           = sim_report.ReservoirTime(rti); % seconds since sim. start
    
    % Get volumes of co2 found in different trapping categories, using
    % heights, topology, trapping structure (ta), etc
    
        % first need to convert saturation into height:
        state       = addCO2HeightData(state, Gt, fluid);
        % state now has .h_free and .h_max fields
        
        % then do some renaming (sol is struct, state is cell)
        sol         = state{:};
        sol.s       = sol.s(:,2);
        sol.h       = sol.h_free;
        fluid.sw    = fluid.res_water;
        fluid.sr    = fluid.res_gas;
        
        % then get normalized values
        %[s, sol.h, sol.h_max]  = normalizeValuesVE(Gt, sol, fluid);

    %ta  = trapAnalysis(Gt, isCellBasedMethod); % done outside function
    vol = volumesVE(Gt, sol, rock2D, fluid, trapstruct);
    
    % total mass injected at the ReservoirTime2plot (or the rti index)
    %totVolInj_m3    = massNow(rti) * 1e9 / rhoCref;
    totVolInj_m3 = sum(vol);                            % TODO: determine why totVolInj_m3 is >> sum(vol) (as if CO2 has exited domain?)

    % since figure is already made, call plotPanelVE() without a varargout
    plotPanelVE(G, Gt, W, sol, t, [vol totVolInj_m3], opts{:});
    
end

hfig = gcf;
hax  = gca;


end

