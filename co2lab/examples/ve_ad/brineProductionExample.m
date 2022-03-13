%% Brine Production Example
% In this example, we use a simple well configuration in the Bjarmeland
% formation (Barents Sea) to demonstrate the use of brine production for
% enhancing storage capacity. Injection and production controls are
% optimized according to an objective function which penalizes leakage and
% pressure-buildup (as done in 'pressureLimitedExample.m').

mrstModule add ad-core ad-props optimization
gravity on;

saveResults = false;


%% Set-up model and perform initial simulation
[model, schedule, initState, seainfo, other] = setupProdExModel( ...
    'itime', 10*year,   'isteps',5,  ...
    'mtime', 100*year,  'msteps',10  );
%%
% Initial simulation:
[init.wellSols, init.states] = simulateScheduleAD(initState, model, schedule, ...
                             'NonLinearSolver', NonLinearSolver('useRelaxation', true));
init.schedule = schedule;

%sameControlTypes = all(strcmpi({schedule.control(1).W.type}, {schedule.control(2).W.type}));


%% Define some parameters
% We use a leakage penalty factor of 5. We set the pressure limit to be 90%
% of the overburden pressure, and the pressure penalty factors will be
% gradually ramped up until the optimal solution is one in which the
% pressure limit is surpassed within a tolerance of 2%.
P_over      = computeOverburdenPressure(model.G, model.rock, seainfo.seafloor_depth, model.fluid.rhoWS);
p_lim_fac   = 0.9;
P_limit     = P_over * p_lim_fac;
cP          = [1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1]; % pressure penalty factor
cL          = 5; % leakage penalty factor


%% Evaluate initial objective (using lowest cP):
% Our objective is to maximize CO2 storage while minimizing long-term
% leakage and pressure buildup. See 'pressureLimitedExample.m' for more
% details on this objective function.
obj_funA = @(wellSols, states, schedule, varargin) ...
            leakPenalizerAtInfinity(model, wellSols, states, ...
                schedule, cL, other.surface_pressure, ...
                model.fluid.rhoWS, other.ta, varargin{:});       
obj_funB = @(states, varargin) ...
            pressurePenalizer(model, states, schedule, cP(1), ...
                P_limit, varargin{:});
           
% We scale the initial objective using obj_funA since obj_funB can be very
% large if the initial guess caused the pressure limit to be greatly
% surpassed.
init.obj_val_steps_A = cell2mat( obj_funA(init.wellSols, init.states, init.schedule) );
init.obj_val_steps_B = cell2mat( obj_funB(init.states) );
init.obj_val_steps = init.obj_val_steps_A - init.obj_val_steps_B;
init.obj_val_total = sum(init.obj_val_steps);
%obj_scaling = abs(init.obj_val_total);
obj_scaling = abs( sum(init.obj_val_steps_A) );

init0 = init; % keep the initial simulation results


%% Define min and max well control values (used to set up box limits)
% Here, we set the limits for the injector wells to be "zero" and some
% multiple of the highest injector rate. The producer wells are limited to
% operate between 50 bars and the initial hydrostatic pressure near the
% well.

% schedule indexes for 'rate' (injector) and 'bhp' (producer) control wells
injWinx = find(strcmpi({schedule.control(1).W.type},'rate'));
prdWinx = find(strcmpi({schedule.control(1).W.type},'bhp'));
prdWcel = [schedule.control(1).W(prdWinx).cells];

% well control limits for injection period
min_wvals( injWinx ) = sqrt(eps);                                       % min well value for injectors, m3/s
min_wvals( prdWinx ) = 50*barsa;                                        % min well value for producers, Pascal
max_wvals( injWinx ) = 3*max([schedule.control(1).W( injWinx ).val]);   % max well value for injectors, m3/s
max_wvals( prdWinx ) = initState.pressure( prdWcel );                   % max well value for producers, Pascals



%% Loop through pressure penalty factors until convergence reached within
% tolerance of 2% of P_limit
sch = schedule;
for r = 1:numel(cP)
    cP_curr = cP(r);

    % Must re-define pressure penalizer with current penalty factor
    obj_funB = @(states, varargin) ...
                pressurePenalizer(model, states, sch, cP_curr, ...
                    P_limit, varargin{:});

    obj_fun = @(wellSols, states, schedule, varargin) ...
                cellSubtract(  obj_funA(wellSols, states, schedule, varargin{:}), ...
                    obj_funB(states, varargin{:})   ); 
                
                
    % Optimize injection rates using the BFGS optimization algorithm.
    % Passing in 'obj_scaling' means an initial simulation will not be
    % executed by optimizeRates. The initial schedule is as per 'sch'.
    [optim, init, history] = optimizeControls(initState, model, sch, min_wvals, max_wvals, ...
                                'obj_fun',                      obj_fun, ...
                                'last_control_is_migration',    true, ...
                                'obj_scaling',                  obj_scaling, ...
                                'lineSearchMaxIt',              5, ...
                                'gradTol',                      1e-3, ...
                                'objChangeTol',                 1e-3);

    
    % Decide whether to ramp up cP or whether convergence reached
    [~, perc_of_Pover_reach] = ...
        findMaxPercentagePlimitReached( optim.states, P_limit, P_over );
    
    if (perc_of_Pover_reach/100 - p_lim_fac) > 0.02
        % use optimized rates as next iteration's initial rates
        sch = optim.schedule;
        % NB: init_scale won't reduce these rates
        % NB: max_wvals won't be in terms of these rates
    elseif r == numel(cP)
        warning('Convergence did not occur before or using last pressure penalty factor.')
    else
        % open system: if p < (plim + tolerance), optimal rates found
        break % exit cp(r) loop
    end
    
end

%% Saving
% add variables to 'other', used for post-processing
other.fluid = model.fluid;
other.initState = initState;
other.leak_penalty = cL;
Gt = model.G;

if saveResults
    dirName = 'brineProductionExampleResults';
    mkdir(dirName);
    save(fullfile(dirName, 'Gt'),       'Gt');
    save(fullfile(dirName, 'optim'),    'optim');
    save(fullfile(dirName, 'init0'),    'init0');
    save(fullfile(dirName, 'init'),     'init');
    save(fullfile(dirName, 'history'),  'history');
    save(fullfile(dirName, 'other'),    'other'); % @@ fluid structure is quite large
end

%%
% <html>
% <p><font size="-1">
% Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.
% </font></p>
% <p><font size="-1">
% This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
% </font></p>
% <p><font size="-1">
% MRST is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% </font></p>
% <p><font size="-1">
% MRST is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% </font></p>
% <p><font size="-1">
% You should have received a copy of the GNU General Public License
% along with MRST.  If not, see
% <a href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses</a>.
% </font></p>
% </html>
