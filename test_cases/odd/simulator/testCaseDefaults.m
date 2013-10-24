function opt = testCaseDefaults


    %% PARAMETERS FOR 50x1 GRID    
   %h_profile_1 = zeros(51,1);  h_profile_1(21:31) = 30;
    temp_grad = 40; % deg/km

    opt =  struct('gravityTheta', (0 * pi / 180), ...
                  'gravityDir',   [1, 0], ...
                  'wellPos',      [25, 1], ... 
                  'h0',           0, ... % h_profile_1
                  'ref_depth',    0, ... % depth for reference press. and temp.
                  'ref_press',    1*atm, ...% pressure at ref. depth
                  'ref_temp',     273.15 + 4 - 0.1 * temp_grad,...
                  'ref_cell',     [1, 1], ... % ref. needed due to inclination
                  'temp_grad',    temp_grad, ...
                  'mu',           [5.36108e-5, 6.5e-4], ...  % [mu_CO2, mu_water]
                  'waterDensity', 1000);


    opt.bcondTypes = {'Pressure', 'Pressure', 'Flux', 'Flux'};    
    
    %% PARAMETERS FOR 20x20 GRID

    % h_profile_2D = zeros(20, 20); h_profile_2D(3:5, 3:5) = 30;    
    % opt =  struct('gravityTheta', (0 * pi / 180), ...
    %               'gravityDir',   [1, 0], ...
    %               'wellPos',      [10 10], ... %[25, 1], ...
    %               'h0',           h_profile_2D(:), ... %0, ...
    %               'p0',           50 * barsa, ...
    %               'p0_pos',       [10, 10], ... %[25, 1], ...
    %               'temperature',  310, ...
    %               'mu',           [5.36108e-5, 6.5e-4]);   % [mu_CO2, mu_water]
    % opt.bcondTypes = {'Pressure', 'Pressure', 'Pressure', 'Pressure'};
    
    
    % options for 'compressible' are: 
    %    - 'incompressible' (incompressible)
    %    - 'horizontal'     (only in xy-direction)
    %    - 'full'           (compressible in x, y and z)

    %% OTHER PARAMETERS
    opt.wellType = {'Rate'};

    % options below are 'incompressible', 'horizontal', 'full'
    opt.compressible = 'full'; % for CO2 (compressible water not yet supported)

    % mass rates, or bottom hole pressures.  One matrix per well.
    opt.schedule = {[1,  0, ...%1e6 * kilo * kilogram / year; ...
                     101, 0]};

    
    opt.simTime   = 100 * year;
    opt.timesteps = 24; 

end


% mu-value for CO2 in Pa/s at 310 K
% mu-value for water in Pa/s at 40 C (313.15 K)

