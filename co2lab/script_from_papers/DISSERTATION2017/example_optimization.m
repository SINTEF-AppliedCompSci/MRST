%% First, run exampleVE up to the point just before simulateScheduleAD starts


% Fixed-rate optimization
min_rate = eps;  % minimum allowed rate should be (practically) zero
max_rate = schedule.control(1).W.val * 2; % maximum allowed rate

schedule.control(2).W.val = min_rate; % no injection during migration

% compute optimized rates, using a leak penalty factor of 10
[optim, init, history] = ...
    optimizeControls(initState, model, schedule, min_rate, max_rate, ...
                  'last_control_is_migration', true, ...
                  'leak_penalty', 10);

% Plotting CO2 staturation for timestep 200 (1100 years after start)
[h, h_max] = upscaledSat2height(optim.states{end}.s(:,2), ...
                                optim.states{end}.sGmax, Gt, ...
                                'pcWG', fluid.pcWG, ...
                                'rhoW', fluid.rhoW, ...
                                'rhoG', fluid.rhoG, ...
                                'p', optim.states{end}.pressure);
plotCellData(Gt.parent, height2Sat(struct('h', h, 'h_max', h_max), Gt, fluid));
colorbar; view(-63, 68); set(gcf, 'position', [531   337   923   356]); axis tight; 


% Plot trapping inventory
ta = trapAnalysis(Gt, false);
reports = makeReports(Gt, {initState, optim.states{:}}, model.rock, model.fluid, ...
                      optim.schedule, [srw, src], ta, []);

h1 = figure; plot(1); ax = get(h1, 'currentaxes');
plotTrappingDistribution(ax, reports, 'legend_location', 'northwest');
set(gcf, 'position', [0 0 1100, 740])
set(gca, 'fontsize', 20);

%% Variable-rate optimization
schedule.control = [repmat(schedule.control(1), 4, 1); ...
                    schedule.control(2)];
schedule.step.control(1:25) = 1;
schedule.step.control(26:50) = 2;
schedule.step.control(51:75) = 3;
schedule.step.control(76:100) = 4;
schedule.step.control(101:200) = 5;

[optim2, init2, history2] = ...
    optimizeControls(initState, model, schedule, min_rate, max_rate, ...
                  'last_control_is_migration', true, ...
                  'leak_penalty', 10);


%% Plotting CO2 staturation for timestep 200 (1100 years after start)
[h, h_max] = upscaledSat2height(optim2.states{end}.s(:,2), ...
                                optim2.states{end}.sGmax, Gt, ...
                                'pcWG', fluid.pcWG, ...
                                'rhoW', fluid.rhoW, ...
                                'rhoG', fluid.rhoG, ...
                                'p', optim2.states{end}.pressure);
plotCellData(Gt.parent, height2Sat(struct('h', h, 'h_max', h_max), Gt, fluid));
colorbar; view(-63, 68); set(gcf, 'position', [531   337   923   356]); axis tight; 

%% Plot trapping inventory
ta = trapAnalysis(Gt, false);
reports = makeReports(Gt, {initState optim2.states{:}}, model.rock, model.fluid, ...
                      optim2.schedule, [srw, src], ta, []);

h1 = figure; plot(1); ax = get(h1, 'currentaxes');
plotTrappingDistribution(ax, reports, 'legend_location', 'northwest');
set(gcf, 'position', [0 0 1100, 740])
set(gca, 'fontsize', 20);
