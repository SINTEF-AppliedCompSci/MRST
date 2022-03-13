%% Simple example of pressure-limited injection for CO2 storage
% This example uses the Bjarmeland formation, which is part of the
% Norwegian Petroleum Directorate's Compiled CO2 Storage Atlas of the
% Norwegian Continental Shelf, found here:
%       http://www.npd.no/en/Publications/Reports/Compiled-CO2-atlas/

% We use a fairly coarsened version of this formation dataset, to reduce
% the required time to run this example. Coarsening is controlled upon call
% to makeBjarmelandModel().

mrstModule add ad-core optimization
saveResults = false;

%% Set up simple example
% The Bjarmeland model is set up with two wells, placed downstream from
% very large structural traps. Their initial rates are set based on the
% trapping capacity of these traps (using prior knowledge of trapping
% structure). However, in this example, we scale down the initial rates
% because it is not possible to operate at such large injection rates
% without comprimising the integrity of the caprock (i.e., the pressure
% limit becomes surpassed).

init_scale = 0.01; % scale down the total injected masses (kg)
qt = init_scale.*[3.7028e12; 2.6198e12];
[model, schedule, initState, seainfo, other] = makeBjarmelandModel( 'qt',qt );

% Initial simulation:
[init.wellSols, init.states] = simulateScheduleAD(initState, model, schedule, ...
                             'NonLinearSolver', NonLinearSolver('useRelaxation', true));
init.schedule = schedule;


%% Set up
% The cell pressure limit is set to be 90% of the overburden pressure,
% which is the weight of the fluid and solid layers lying on top of the
% formation we are injecting into. Cell pressure which surpasses this limit
% will be penalized by the objective function, using an optimization
% strategy in which the penalty factor is gradually ramped up as the
% iterations of the optimization scheme progress. Convergence to the
% optimal solution is reached once cell pressure surpasses the pressure
% limit by a tolerance of 2%. Leakage is also penalized by the objective
% function, and here we use a penalty factor of 5.

P_over = computeOverburdenPressure(model.G, model.rock, seainfo.seafloor_depth, model.fluid.rhoWS);
p_lim_fac   = 0.9;
P_limit     = P_over * p_lim_fac;
cP = [1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1]; % pressure penalty factor
cL = 5; % leakage penalty factor
rate_lim_fac = 2; % max well controls are defined in terms of non-scaled qt


%% Evaluate initial objective (using lowest cP):
%   Basic form: J = Ja - Jb, where
%       Ja = MI - cL*ML
%       Jb = cP*sum( (max(0,P-P_limit)^2 )
% Here, we set up the objective function, which is constructed to penalize
% leakage and pressure buildup. We compute the objective function value for
% the initial simulation. This value will be used to scale the objective
% values computed within the main optimization routine.

obj_funA = @(wellSols, states, schedule, varargin) ...
            leakPenalizerAtInfinity(model, wellSols, states, ...
                schedule, cL, other.surface_pressure, ...
                model.fluid.rhoWS, other.ta, varargin{:});       
obj_funB = @(states, varargin) ...
            pressurePenalizer(model, states, schedule, cP(1), ...
                P_limit, varargin{:});
           
% Scale initial objective using obj_funA since obj_funB can be very large
init.obj_val_steps_A = cell2mat( obj_funA(init.wellSols, init.states, init.schedule) );
init.obj_val_steps_B = cell2mat( obj_funB(init.states) );
init.obj_val_steps = init.obj_val_steps_A - init.obj_val_steps_B;
init.obj_val_total = sum(init.obj_val_steps);
%obj_scaling = abs(init.obj_val_total);
obj_scaling = abs( sum(init.obj_val_steps_A) );

init0 = init; % keep the initial simulation results


%% More set up
% We set how small and how high the well rates are allowed to be, which
% constrains our optimization problem within a "box".

% Define limits and scaling based on initial simulation:
min_wvals       = sqrt(eps) * ones(numel(schedule.control(1).W), 1);
max_wvals       = rate_lim_fac * max( [schedule.control(1).W.val]' ./ init_scale ) ...
                    * ones(numel(schedule.control(1).W), 1);
                

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
                
                
    % Optimize injection rates
    % using the BFGS optimization algorithm. Passing in 'obj_scaling' means
    % an initial simulation will not be executed by optimizeRates. The
    % initial schedule is as per 'sch'.
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
    dirName = 'pressureLimitedExampleResults';
    mkdir(dirName);
    save(fullfile(dirName, 'Gt'),       'Gt');
    save(fullfile(dirName, 'optim'),    'optim');
    save(fullfile(dirName, 'init0'),    'init0');
    save(fullfile(dirName, 'init'),     'init');
    save(fullfile(dirName, 'history'),  'history');
    save(fullfile(dirName, 'other'),    'other'); % @@ fluid structure is quite large
end

%% Post-processing
% The following figures will be generated:
%   1. well placement with trapping structure
%   2. initial co2 saturation
%   3. optimized co2 saturation
%   4. overburden pressure
%   5. max fraction of overburden pressure reached
%   6. initial vs optimized well rates
%   7. initial vs optimized trapping inventory (with forecast curves)
%   101: breakdown of forecast given initial rates
%   102: breakdown of forecast given optimized rates
close all
postProcessExample(Gt, init0, optim, other);

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
