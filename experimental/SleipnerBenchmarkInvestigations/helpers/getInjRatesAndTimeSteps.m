function [ schedule, var ] = getInjRatesAndTimeSteps( varargin )


    % Default rate is IEAGHG
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


end


