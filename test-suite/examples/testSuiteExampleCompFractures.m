%% Test suite example: compositional flow in fractured reservoir
% In this example, we consider a 2D example with compositional flow in a
% fractured reservoir. We show how to set up and simulate the problem with
% a fully implicit and sequential implict solver, and compare the results

%% Add required modules
mrstModule add test-suite
mrstModule add ad-core ad-props ad-blackoil
mrstModule add compositional
mrstModule add linearsolvers sequential
mrstModule add mrst-gui
mrstVerbose on
checkHashSettings();

%% Load test case
testFI = TestCase('fractures_compositional');

%% Plot test case
lperm = log10(convertTo(testFI.model.rock.perm(:,1), milli*darcy));
testFI.plot(lperm);
plotGrid(testFI.model.G, 'faceColor', 'none', 'edgeAlpha', 0.2);
[hc, hh] = mrstColorbar(10.^lperm, 'south', true);
colormap(pink);

%% Simulate FI problem
problemFI = testFI.getPackedSimulationProblem('useHash', true);
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
tsol = AMGCLSolverAD('tolerance'     , 1e-4        , ...
                     'preconditioner', 'relaxation', ...
                     'relaxation'    , 'ilu0'      , ...
                     'block_size'    , ncomp       );
ordering = getCellMajorReordering(model.G.cells.num, ncomp);
tsol.equationOrdering = ordering;
tsol.variableOrdering = ordering;
% Set solvers
testSI.model.pressureNonLinearSolver.LinearSolver  = psol;
testSI.model.transportNonLinearSolver.LinearSolver = tsol;

%% Simulate SI problem
problemSI = testSI.getPackedSimulationProblem('useHash', true);
simulatePackedProblem(problemSI);

%% Get results and compare FI to SI
[wellSolsFI, statesFI, reportsFI] = getPackedSimulatorOutput(problemFI);
[wellSolsSI, statesSI, reportsSI] = getPackedSimulatorOutput(problemSI);
% To compare the results, we use the function `compareStructs`. This 
compare = @(sFI, sSI) compareStructs(sFI, sSI, 'relative'      , true ,  ...
                                               'includeStructs', false, ...
                                               'verbose'       , false);
statesDiff = cellfun(@(sFI, sSI) compare(sFI, sSI), statesFI, statesSI, ...
                                                'UniformOutput', false);
                                            
%% Plot results
testFI.plot(statesFI  , 'Name', 'Fully implicit'     ); colormap(bone);
testSI.plot(statesSI  , 'Name', 'Sequential implicit'); colormap(bone);
testFI.plot(statesDiff, 'Name', 'Difference'         ); colormap(flipud(bone));
