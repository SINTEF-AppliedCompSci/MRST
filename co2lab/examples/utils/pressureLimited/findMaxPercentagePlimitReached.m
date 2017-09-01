function [perc_of_plim_reach, perc_of_Pover_reach, cinxMax, cinxViolated] = ...
    findMaxPercentagePlimitReached( states, plim, P_over, varargin )
% Determine maximum percentage of pressure limit (or overburden pressure)
% that was reached at some time and location during the simulation


% states    - cell array containing pressure fields
% plim      - pressure limit (some percentage of P_over)
% P_over    - overburden pressure
% cinxMax   - cell index where maximum percentage of plim reached

    opt.outputOn = true;
    opt = merge_options( opt, varargin{:} );

    % Maximum pressure encountered during simulation
    tmp = [];
    for j=1:numel(states)
        tmp = [tmp, states{j}.pressure];
    end
    %maxP_encountered = max(tmp')'; % cell array, 1:Gt.cells.num
    maxP_encountered = max(tmp,[],2); % takes max value along rows instead 
                                      % of column (to handle possibility
                                      % that states is a single column
                                      % array, not a cell array of many
                                      % column arrays)
    
    % Fraction of plim reached.
    frac = maxP_encountered./plim;

    % Find the largest fraction of plim that was reached. If largest
    % fraction is 1 or above, plim has been reached or surpassed.
    % Otherwise, plim has not been reached.
    [~, cinxMax] = max(frac);
    
    cinxViolated = find(frac > 1);
    
    if opt.outputOn
        fprintf('Pressure reached %4.3f percent of plim, in cell %d, ... \n', ...
            maxP_encountered(cinxMax)/(plim(cinxMax))*100, cinxMax)
        fprintf('...which is %4.3f percent of the overburden pressure.\n', ...
            maxP_encountered(cinxMax)/P_over(cinxMax)*100)
    end
    perc_of_plim_reach = maxP_encountered(cinxMax)/(plim(cinxMax))*100;
    perc_of_Pover_reach = maxP_encountered(cinxMax)/P_over(cinxMax)*100;




end

