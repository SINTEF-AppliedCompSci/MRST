%% Resolution of trailing waves with high-resolution schemes
% The example runs a simple displacement for which the c-wave should be a
% shock (Todd-Longstaff mixing parameter w=0.85). The purpose of the
% example is to demonstrate that using a high-resolution spatial
% discretization significantly improves the resolution of the trailing
% chemical wave.
mrstModule add ad-core ad-blackoil ad-eor ad-props deckformat

%% Run fine-grid simulation
gravity reset on;
bookdir = getDatasetPath('eor_book_ii');
fn = fullfile(bookdir,'trailingWaves','MixPar400.DATA');
[state0, model, schedule, nlsolver] = initEclipseProblemAD(fn);
model.fluid.mixPar = 0.85;
[~, states] = simulateScheduleAD(state0, model, schedule, 'NonLinearSolver', nlsolver);

clf
subplot(1,2,1), hold all
plot(model.G.cells.centroids(:,1),states{end}.s(:,1),'LineWidth',2);
subplot(1,2,2), hold all
plot(model.G.cells.centroids(:,1),states{end}.cp(:,1),'LineWidth',2);

%% Run simulation with 50 cells
fn = fullfile(bookdir,'trailingWaves','MixPar50.DATA');
[state0, model, schedule, nlsolver] = initEclipseProblemAD(fn);
model.fluid.mixPar = 0.85;
[~, states] = simulateScheduleAD(state0, model, schedule, 'NonLinearSolver', nlsolver);

pargs = {'-','Marker','.','MarkerSize', 14, 'LineWidth',1};
subplot(1,2,1)
plot(model.G.cells.centroids(:,1),states{end}.s(:,1),pargs{:});
subplot(1,2,2)
plot(model.G.cells.centroids(:,1),states{end}.cp(:,1),pargs{:});

%% Run the same with WENO and with WENO + AIM
model = setWENODiscretization(model);
[~, states] = simulateScheduleAD(state0, model, schedule, 'NonLinearSolver', nlsolver);
subplot(1,2,1)
plot(model.G.cells.centroids(:,1),states{end}.s(:,1),pargs{:});
subplot(1,2,2)
plot(model.G.cells.centroids(:,1),states{end}.cp(:,1),pargs{:});

model = setTimeDiscretization(model, 'aim', 'saturationCFL', 0.75);
[~, states] = simulateScheduleAD(state0, model, schedule, 'NonLinearSolver', nlsolver);
subplot(1,2,1)
plot(model.G.cells.centroids(:,1),states{end}.s(:,1),pargs{:});
subplot(1,2,2)
plot(model.G.cells.centroids(:,1),states{end}.cp(:,1));
legend('SPU: \Delta x = 1 m','SPU: \Delta x = 8 m', ...
    'WENO: \Delta x = 8 m', 'WENO+AIM: \Delta x = 8 m');