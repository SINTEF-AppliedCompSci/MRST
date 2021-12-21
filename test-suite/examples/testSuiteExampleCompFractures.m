%% Add required modules
mrstModule add test-suite
mrstModule add ad-core ad-props ad-blackoil
mrstModule add compositional
mrstModule add sequential
mrstModule add mrst-gui
mrstVerbose on

%% Load test case
testFI = TestCase('fractures_compositional');

%% Plot test case
lperm = log10(testFI.model.rock.perm(:,1)/(milli*darcy));
testFI.plot(lperm);
plotGrid(testFI.model.G, 'faceColor', 'none', 'edgeAlpha', 0.2);
[hc, hh] = mrstColorbar(10.^lperm, 'south', true);
colormap(pink);

%% Simulate FI problem
problemFI = testFI.getPackedSimulationProblem();
simulatePackedProblem(problemFI);

%% Set sequential implicit formulation (SI)
testSI = testFI;
testSI.model = getSequentialModelFromFI(testSI.model);
% For the pressure subproblem, we use an algebraic multigrid solver with
% incomplete LU factorization with zero fill-in (ILU0)
psol = AMGCLSolverAD('tolerance', 1e-4);
% For the transport subproblem, we use a Krylov solver with ILU(0)
% preconditioning
% Get number of components
model = testFI.model.validateModel();
ncomp = model.getNumberOfComponents();
tsol  = AMGCL_CPRSolverAD('block_size'    , ncomp       ,...
                          'preconditioner', 'relaxation', ...
                          'relaxation'    , 'ilu0'      , ...
                          'tolerance'     , 1e-4        );
% Set solvers
testSI.model.pressureNonLinearSolver.LinearSolver  = psol;
testSI.model.transportNonLinearSolver.LinearSolver = tsol;

%% Simulate SI problem
problemSI = testSI.getPackedSimulationProblem();
simulatePackedProblem(problemSI);

%% Get results and compare FI to SI
[wellSolsFI, statesFI, reportsFI] = getPackedSimulatorOutput(problemFI);
[wellSolsSI, statesSI, reportsSI] = getPackedSimulatorOutput(problemSI);
% To compare the results, we use the function `compareStructs`. This 
compare = @(sFI, sSI) compareStructs(sFI, sSI, 'relative'      , true ,  ...
                                               'includeStructs', false);
statesDiff = cellfun(@(sFI, sSI) compare(sFI, sSI), statesFI, statesSI, ...
                                                'UniformOutput', false);
                                            
%% Plot results
test.plot(statesFI  , 'Name', 'Fully implicit'     ); colormap(bone);
test.plot(statesSI  , 'Name', 'Sequential implicit'); colormap(bone);
test.plot(statesDiff, 'Name', 'Difference'         ); colormap(bone);