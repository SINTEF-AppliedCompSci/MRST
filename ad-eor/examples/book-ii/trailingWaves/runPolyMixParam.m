%% Resolution of trailing waves
% The example runs a simple displacement for various mixing parameters in
% the Todd-Longstaff method. In all cases, the c-wave should be a
% discontinuity (contact discontinuity for w=1 and shock otherwise).
% The purpose of the example is to demonstrate that the resolution of the
% c-wave depends on the degree of shelf sharpening in the wave, which
% increases as w decreases.
%
% Altogether: there are four different input files you could use:
%  - MixPar25.DATA:     25 grid cells, dx=16 m
%  - MixPar50.DATA:     50 grid cells, dx= 8 m
%  - MixPar400.DATA:   400 grid cells, dx= 1 m
%  - MixParAds25.DATA:  25 grid cells, with adsorption and RRF
mrstModule add ad-core ad-blackoil ad-eor ad-props deckformat

%% Load and setup simulation framework
gravity reset on;
bookdir = getDatasetPath('eor_book_ii');
fn = fullfile(bookdir,'trailingWaves','MixPar25.DATA');
[state0, model, schedule, nlsolver] = initEclipseProblemAD(fn);

%% Run simulations with various mixture parameters
mixPar = [1:-.025:.9 .85 .8 .6 0];
col = jet(numel(mixPar));
set(gcf,'Position',[370 530 1040 224]);
for i=1:numel(mixPar)
    model.fluid.mixPar = mixPar(i);
    [~, states] = simulateScheduleAD(state0, model, schedule, ...
        'NonLinearSolver', nlsolver);
    
    subplot(1,2,1), hold all
    plot(model.G.cells.centroids(:,1),states{end}.s(:,1), ...
        'LineWidth',.5+1.5*ismember(mixPar(i),[0 .9 1]),'Color',col(i,:));
    subplot(1,2,2), hold all
    plot(model.G.cells.centroids(:,1),states{end}.cp(:,1), ...
        'LineWidth',.5+1.5*ismember(mixPar(i),[0 .9 1]),'Color',col(i,:));
    drawnow
end
legend(cellfun(@(x) sprintf('\\omega=%.3f', x), num2cell(mixPar),...
    'UniformOutput', false));
subplot(1,2,1), set(gca,'FontSize',12), axis([0 400 0.1 0.85]);
subplot(1,2,2), set(gca,'FontSize',12), axis([0 400 -0.1 3.1]);
set(gca,'Position',get(gca,'Position')-[.05 0 0 0]);