function [ schedule, var ] = getInjRatesAndTimeSteps( varargin )


    % Default rates is IEAGHG
    opt.ratecase = 'SPE';
    
    % for limiting the number of injection years
    opt.num_years   = [];
    
    

    opt = merge_options(opt, varargin{:});
    
    
    %% Select well injection rate data:
    fprintf(['Your rate case is set to ' opt.ratecase '.\n'])
    switch opt.ratecase 
        % Note, inj_rates are in terms of reservoir rates (i.e., the
        % volumetric rate of CO2 entering layer 9, not the volumetric
        % surface rate). Seismic imaging provided estimates of how much CO2
        % accumlated in the pore space of layer 9. These volumes were
        % likely converted into a mass using an infered CO2 density, and
        % then into a surface rate using the CO2 density at the surface.
        % Specifying the inj_rates in terms of reservoir volume instead of
        % surface volume allows one to test other CO2 densities without the
        % need to modify a surface volume injection rate.
        
        case 'SPE'
            % See Singh et al 2010 for more info about how they determined
            % these rates. Note: the injection rates were reported as
            % surface rates. Both volume and mass were given, thus surface
            % density can be calculated (=1.87 kg/m3). The CO2 density at
            % reservoir conditions was reported as 760 kg/m3.
            var.inj_year   = [1999; 2000; 2001; 2002; 2003; 2004; 2005; 2006; 2007; 2008; 2009];
            var.inj_rates  = [2.91e4; 5.92e4; 6.35e4; 8.0e4; 1.09e5; 1.5e5; 2.03e5; 2.69e5; 3.47e5; 4.37e5; 5.4e5] .* meter^3/day;
            % inj_rates is in meter^3/s
            % Convert to rate at reservoir conditions
            var.inj_rates  = var.inj_rates.*(1.87/760);       
        case 'original' 
            % See "Injection rates Layer 9.xls" under
            % co2lab/data/sleipner/original for more info about these
            % rates. Note: the CO2 density at reservoir conditions was
            % reported as 695.89 kg/m3, and the surface density 1.87 kg/m3.
            [ var.inj_year, var.inj_rates ] = getSleipnerOriginalInjectionRates();
            % inj_rates is in meter^3/s   
        otherwise
            error('The injection rate option was either invalid or not selected.')     
    end
    
    
    % limit the number of injection years (optional)
    if ~isempty(opt.num_years)
        fprintf('\n Only using first %d years of injection data. \n', opt.num_years)
        var.inj_year  = var.inj_year(1:opt.num_years);
        var.inj_rates = var.inj_rates(1:opt.num_years);
    end
    
    
    
    
    %% Compute time step size for injection period
    
    % each inj time is 1 year (converted into seconds), since yearly rates
    % have been given
    var.inj_time    = ( ones( numel(var.inj_year), 1 ) ) * year;
    
    % default step size is 1 step per year
    var.inj_steps   = ones( numel(var.inj_year), 1 );
    
    % however, specify more steps to use during the first several years:
    var.inj_steps(1:10,1)   = [ 8; 4; 4; 2; 2; 2; 2; 2; 2; 2 ];
    if ~isempty(opt.num_years)
        var.inj_steps  = var.inj_steps(1:opt.num_years);
    end
    
    % time step sizes to use in each year:
    var.dTi     = var.inj_time ./ var.inj_steps;

    
    % Put time step sizes into schedule ---
    % time step size value:
    schedule.step.val = [];
    for i = 1:numel(var.inj_steps)
        schedule.step.val = [schedule.step.val; ...
            repmat( ones(var.inj_steps(i), 1) * var.dTi(i) , [numel(var.inj_rates(i)) 1] )]; 
    end
    %schedule.step.val = [schedule.step.val(1); schedule.step.val(1:end)];

    % time step control (i.e., will be the corresponding control well):
    schedule.step.control   = [];
    for i = 1:numel(var.inj_rates)
        schedule.step.control = [schedule.step.control; ...
            ones(var.inj_steps(i), 1) * i ];
    end
    

    
    
%     % TODO: implement for greater flexibility
%     timeStepIntervals = [1999:0.25:2000, 2000.5:0.5:2002.5, 2003.5:1:2009.5]'; %, 2010]';
%     for i = 1:numel(timeStepIntervals)-1
%         dTi(i,1) = timeStepIntervals(i+1) - timeStepIntervals(i);
%     end
%     dTi = dTi * year; % timestep size in seconds
%     
%     % Inj rates corresponding to dTi
%     inj_rates_per_dTi       = zeros(numel(var.inj_rates),1);
%     inj_rates_per_dTi(1)    = var.inj_rates(1);
%     for i = 2:numel(dTi)
%         inj_rates_per_dTi(i,1) = (var.inj_rates(i-1) + var.inj_rates(i))/2;
%     end
%     %inj_rates_per_dTi(end+1,1) = var.inj_rates(end);
%     
%     plot(timeStepIntervals(1:end-1), inj_rates_per_dTi, 'x')
%     hold on;
%     plot(var.inj_year, var.inj_rates, 'o')
%     grid; legend('inj rate interpolated','inj rate data', 'Location','northwest')
    
    
%     %% Add well rates to Schedule
%     % injection period only
%     for i = 1:numel(var.inj_rates)
%         schedule.control(i).W = addWell([], Gt.parent, rock2D, ...
%             wellCellIndex, 'name', sprintf('W%i', i), 'Type', 'rate', ...
%             'Val', var.inj_rates(i), 'comp_i', [0 1]);
%             % inj_rate should be mass rate / fluid.rhoGS
%     end
% 
% 
%     %% Add time steps to Schedule
%     % TODO: implement dynamic time-stepping
%     %istepvec                = repmat( ones(inj_steps, 1) * dTi , [numel(var.inj_rates) 1] );   
%     schedule.step.val       = dTi; %istepvec;
%     schedule.step.control   = [];
%     for i = 1:numel(schedule.control)
%         % injection period only
%         schedule.step.control = [schedule.step.control; ones(inj_steps, 1) * i];
%     end


end


