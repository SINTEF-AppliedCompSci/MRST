function [perc_of_plim_reach, perc_of_Pover_reach, cinx] = ...
    report_maxPercentage_plim_reached( states, plim, P_over, varargin )
% Determine maximum percentage of pressure limit that was reached at
% some time and location during the simulation

% P_over * plim_factor = plim
% states.pressure

    opt.outputOn = true;
    opt = merge_options( opt, varargin{:} );

    % Maximum pressure encountered during simulation
    tmp = [];
    for j=1:numel(states)
        tmp = [tmp, states{j}.pressure];
    end
    maxP_encountered = max(tmp')'; % cell array, 1:Gt.cells.num
    
    % Difference from pressure limit.
    % Negative values mean: plim > maxP_encountered
    % Positive values mean: maxP_encountered > plim
    diff = maxP_encountered - plim;
    
    % If any values are positive (i.e., 0 and above), plim has been reached
    % or surpassed. Thus find the largest surpassing value (i.e., the
    % largest positive value)
    if any(diff >= 0)
        [~, cinx] = max(diff);
    end
    % If no positive values, plim has not been surpassed. Find how close
    % the pressure came to plim by taking the minimum of the abs(diff).
    if all(diff < 0)
        [~, cinx] = min(abs(diff));
    end
    
    if opt.outputOn
        fprintf('Pressure reached %4.3f percent of plim, in cell %d, ... \n', ...
            maxP_encountered(cinx)/(plim(cinx))*100, cinx)
        fprintf('...which is %4.3f percent of the overburden pressure.\n', ...
            maxP_encountered(cinx)/P_over(cinx)*100)
    end
    perc_of_plim_reach = maxP_encountered(cinx)/(plim(cinx))*100;
    perc_of_Pover_reach = maxP_encountered(cinx)/P_over(cinx)*100;




end

