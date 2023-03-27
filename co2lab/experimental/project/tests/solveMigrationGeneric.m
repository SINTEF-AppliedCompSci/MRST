function [states, timesteps] = solveMigrationGeneric(Gt, sol, W, bc, fluid, solver, varargin)
opt = struct('T_injection',  100*year,  ...
            'T_migration',  1000*year, ...
            'dt_min', 0.01*year,...
            'dt_max', 10*year,...
            'dt',       [], ...
            'dynamicTimestepping', true, ...
            'cutTimesteps', false, ...
            'dtplot',   1*year, ...
            'dtPost_min', 1*year,... 
            'dtPost_max', 25*year,...
            'doPlot',   true, ...
            'plotter',  @(sol, t) basicPlotter(Gt, sol, t), ...
            'Verbose', mrstVerbose...
            );

opt = merge_options(opt, varargin{:});


Gt.grav_pressure = @(g, omega) gravPressureVE_s(g, omega);
Gt.primitives = @primitivesMimeticVE_s;

gravity on

T          = opt.T_injection + opt.T_migration;
stopInject = opt.T_injection;
if isempty(opt.dt)
    dt         = opt.dt_max;
else
    dt         = opt.dt;
end
dt_min     = opt.dt_min;
dt_max     = opt.dt_max;
dtpost_min = opt.dtPost_min;
dtpost_max = opt.dtPost_max;
dtplot     = opt.dtplot;



t = 0;
timeSincePlot = 0;

dispif(1,'\nSimulating %d years of injection', convertTo(stopInject,year));
dispif(1,' and %d years of migration\n', convertTo(T-stopInject,year));
dispif(1,'Time: %4d years', convertTo(t,year));
stime=tic;

moduleCheck('ad-fi');

isMigrating = false;

dt_history = [];
timesteps = [];
states = [];
its = 0;
while t<T
    switchedToMigration = false;
    sol_previous = sol;
    
    if opt.dynamicTimestepping
        [dt, dt_history] = simpleStepSelector(dt_history, dt, its,...
                                    'dt_min', dt_min, 'dt_max', dt_max);
    end
    
    if (t + dt) > opt.T_migration && ~isMigrating
        dt = opt.T_migration - t;
        
        isMigrating = true;
        switchedToMigration = true;
        dt_max = dtpost_max;
        dt_min = dtpost_min;
        dt_history = [];
    end
    
    [sol, its, conv] = solver(sol, W, bc, dt);
    
    if ~conv.converged && opt.cutTimesteps
        dt = dt/2;
        sol = sol_previous;
        if dt <= dt_min
            warning('Cutting of timestep violated dt_min -> Impossible to converge');
        else
            continue
        end
    end
    timesteps = [timesteps; dt];  %#ok
    states    = [states;    sol]; %#ok
    % Plotting stuff
    timeSincePlot = timeSincePlot + dt;
    if timeSincePlot > dtplot
        % Do plot
        opt.plotter(sol, t + dt);
        timeSincePlot = 0;
    end
    t = t + dt;
    dispif(1,'Time: %s\n', formatTimeRange(t));
    
    if switchedToMigration
        % disable wells
        W = [];
        its = 0;
    end
end

cputime=toc(stime);%#ok

end

function basicPlotter(G, sol, t)
    persistent fig h;
    if isempty(fig) || ~ishandle(fig)
        fig = figure;
    end
    set(0, 'CurrentFigure', fig);
    
    if ~isempty(h) && ishandle(h) 
        delete(h);
    else
        plotGrid(G, 'facecolor', 'none');
    end
    
    % h = plotCellData(G, sol.s, sol.s > sqrt(eps));
    % caxis([0 1])
    h = plotCellData(G, sol.h, sol.h > sqrt(eps));
    plotGrid(G, 'facecolor', 'none', 'edgealpha', .01)
    fastRotateButton
    title(formatTimeRange(t));
    
    colorbar
    drawnow
end