
function system = prepareSimulationSystem(tcase)

    %% Setting up 's' helper structure
    % calling `setupSimComp` instead of `setupSimCompVe`, since there's at least
    % one bug in the latter, and since it's apparently prepared for dealing with
    % saturations/pore volume rather than porosity.  Instead, I call the regular
    % `setupSimComp function, but modify the permeability of the rock first, in
    % order to make it not the vertical average, but the vertical sum.
    rock = tcase.rock;
    rock.perm = rock.perm .* tcase.Gt.cells.H;
    s = setupSimComp(tcase.Gt, rock);
    s.top_temp = tcase.top_temp;
    s.temp_grad = tcase.temp_grad;
    
    %% Setting up system
    system.s = s;                         % storing helper structure in system

    system.nonlinear.maxIterations = 25;  % queried inside 'stepCO2'
    system.nonlinear.changeWells = false; % queried in 'solvefiADI'
    system.nonlinear.bhpcontrols = false; % queried in 'solvefiADI'
    system.nonlinear.relaxation  = 1;     % queried in 'solvefiADI' but not used!
    system.nonlinear.tolMB = 1e-3; %1e-3;%1e-7;        % queried in 'getCO2Convergence'
    system.nonlinear.tolCNV = 1e-2; %1e-2;%1e-3;       % queried in 'getCO2Convergence'

    system.stepFunction = ...
       @(x0, x, meta, dt, W, G, sys) stepCO2(x0, x, meta, dt, W, G, sys, tcase.fluid);

    system.getEquations = @(state0, state, dt, G, W, s, f) ...
        eqsfiCO2water(state0, state, dt, G, W, s, f, tcase.bc, tcase.gravity, ...
                      tcase.top_temp, tcase.temp_grad);
end
