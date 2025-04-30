clear;

mrstModule add ad-core ad-props compositional ecpa dual-porosity dual-porosity-permeability deckformat
gravity on

%% Choose between the dual porosity or dual porosity-dual permeability models
model_name = 'DP';
% model_name = 'DPDP';

%% Set up grid
G = cartGrid([1, 1, 1], [1, 1, 1]*meter);
G = computeGeometry(G);

%% Set up rock properties
rock_matrix = makeRock(G, 0.5*darcy, 0.2);
rock_fracture = makeRock(G, 1*darcy, 0.01);

%% Set up fluid
mixture = ECPATableCompositionalMixture({'Water','CarbonDioxide'});
T = 273.15 + 30; p = 5E6;  z = [1, 0];     % p, T, z
bic = eCPAreadBinaryInteraction(mixture,T(1));
ap = eCPAreadAssociationParameter(mixture, T(1));
mixture = setBinaryInteraction(mixture, bic);
mixture = setAssociationParameter(mixture, ap);

eCPA = ECPAEquationOfStateModel([], mixture, 'eCPA');
[~, ~, ~,~, ~,rhoOSC,rhoGSC] = eCPAstandaloneFlash(p, T, [0 1], eCPA);
fluid_matrix = initSimpleADIFluid('phases', 'OG', 'blackoil', false, 'rho', [1000,50]);
fluid_fracture = initSimpleADIFluid('phases', 'OG', 'blackoil', false, 'rho', [1000,50]);

%% Set the DP model. Here, a single-phase model is used. Rock and fluid
% are sent in the constructor as a cell array {fracture,matrix}
diagonal_backend = DiagonalAutoDiffBackend('modifyOperators', true, 'rowMajor', true);
if strcmp(model_name, 'DP')
    % Initialize the 2 phase dual porosity model using the same fluid model for
    % the matrix and fracture regions
    natural= ECPANaturalVariablesCompositionalModelDP(G, {rock_fracture,rock_matrix},...
                                                        {fluid_fracture,fluid_matrix},...
                                                        mixture,...
                                                        'water',false, 'oil',true, 'gas', true,...
                                                        'liquidPhase', 'O', 'vaporPhase', 'G',...
                                                        'AutoDiffBackend', diagonal_backend);
    overall = ECPAOverallCompositionCompositionalModelDP(G, {rock_fracture,rock_matrix},...
                                                            {fluid_fracture,fluid_matrix},...
                                                            mixture,...
                                                            'water',false, 'oil',true, 'gas', true,...
                                                            'liquidPhase', 'O', 'vaporPhase', 'G',...
                                                            'AutoDiffBackend', diagonal_backend);
elseif strcmp(model_name, 'DPDP')
    % Initialize the 2 phase dual porosity-dual permeability model using the same fluid model for
    % the matrix and fracture regions
    natural= ECPANaturalVariablesCompositionalModelDPDP(G, {rock_fracture,rock_matrix},...
                                                            {fluid_fracture,fluid_matrix},...
                                                            mixture,...
                                                            'water',false, 'oil',true, 'gas', true,...
                                                            'liquidPhase', 'O', 'vaporPhase', 'G',...
                                                            'AutoDiffBackend', diagonal_backend);
    overall = ECPAOverallCompositionCompositionalModelDPDP(G, {rock_fracture,rock_matrix},...
                                                            {fluid_fracture,fluid_matrix},...
                                                            mixture,...
                                                            'water',false, 'oil',true, 'gas', true,...
                                                            'liquidPhase', 'O', 'vaporPhase', 'G',...
                                                            'AutoDiffBackend', diagonal_backend);
end
overall = overall.validateModel();
natural = natural.validateModel();
%% Setting transfer function. This step is important to ensure that fluids
% will be transferred from fracture to matrix (and vice-versa). There are 
% several options of shape factor (see folder
% transfer_models/shape_factor_models/) that could be set in this
% constructor. The second argument is the matrix block size. Another
% possible transfer function to be used in this model would be:
%       model.transfer_model_object = EclipseTransferFunction();
overall.transfer_model_object = KazemiSinglePhaseTransferFunction('KazemiShapeFactor',...
                                                                 [1,1,1]);
natural.transfer_model_object = KazemiSinglePhaseTransferFunction('KazemiShapeFactor',...
                                                                 [1,1,1]);
%% Initializing state
state = eCPAinitCompositionalState(overall, p, T, [], z, true);
state = overall.validateState(state);
initGuess = state;
initGuess = overall.computeFlash(initGuess);
%% Boundary conditions
bc = fluxside([], G, 'xmin', 1/day, 'sat', [0, 1]);    % Flux
bc = pside(bc, G, 'xmax', p, 'sat', [0, 1]);           % Standard bc
bc.components = repmat([0, 1], numel(bc.face), 1);     % Boundary z
bc.value_matrix = bc.value;
%% Solver
nls = NonLinearSolver('reportLevel', 3,'maxIterations', 100);

dt = 100*day;
[solNat1, reportsNat1] = nls.solveTimestep(state, dt, natural, ...
                            'bc', bc, 'initialGuess', initGuess);
[solMole1, reportsMole1] = nls.solveTimestep(state, dt, overall, ...
                            'bc', bc, 'initialGuess', initGuess);

%% Extract data from the output states
getStates = @(reports) cellfun(@(x) x.state, reports.StepReports{1}.NonlinearReport, 'UniformOutput', false);
% Get the outputs
getZ = @(states) [initGuess.components(1); cellfun(@(x) x.components(:, 1), states)];
getS = @(states) [initGuess.s(:, 1); cellfun(@(x) x.s(:, 1), states)];
getP = @(states) [initGuess.pressure(1); cellfun(@(x) x.pressure, states)] - p;
natStates = getStates(reportsNat1);
moleStates = getStates(reportsMole1);

zn = getZ(natStates);
sn = getS(natStates);
dpn = getP(natStates);

zm = getZ(moleStates);
sm = getS(moleStates);
dpm = getP(moleStates);

%% Plot the convergence histories for both models
% We make two different plots of the convergence: in (z_w,dp) space and in
% (S_l,dp) space. These plots are imposed on top of a contour map of S.
ns = 50;
ps = linspace(min(min(dpn), min(dpm)), 1.1*max(max(dpn), max(dpm)), ns);
zs = linspace(0.01, 1, ns);

[DP, Z] = meshgrid(ps, zs);
zz = reshape(Z, [], 1);

[L, ~, ~, Z_L, Z_V] = standaloneFlash(DP(:) + p, T, [zz, 1 - zz], overall.EOSModel);
Lg = reshape(L, ns, ns);
S  = reshape(L.*Z_L./(L.*Z_L + (1-L).*Z_V), ns, ns);
maxChange =0.1;
for i = 1:2
    if i == 1
        xn = zn; xm = zm; XX = Z; l = 'z_{water}';
    else
        xn = sn; xm = sm; XX = S; l = 'S_L';       
    end
    figure(i); clf; hold on
    contourf(XX, DP, S, 10);
    c = lines(4);
    c1 = c(2, :);
    c2 = c(3, :);
    if maxChange < 0.1
        style = '-';
    else
        style = '-o';
    end
    h1 = plot(xn, dpn, style, 'color', c1, 'linewidth', 1.2, ...
        'MarkerFaceColor', 'r', 'MarkerSize', 8); hold on
    h2 = plot(xm, dpm, style, 'color', c2, 'linewidth', 1.2, ...
        'MarkerFaceColor', 'y');
    h3 = plot(xm(end), dpm(end), 'x', 'color', [1, 1, 1]*0.7, ...
        'MarkerSize', 10, 'linewidth', 1.2);
    legend([h1; h2; h3], 'Natural variables', 'Overall composition', ...
        'Solution', 'location', 'southeast')
    ylabel('\Delta p [Pa]')
    xlabel(l)

    cmap = interp1([1; 100], [.1, .1, 1; 0.9, 0.9, 1], linspace(1, 100, 100));
    colormap(flipud(cmap))
end