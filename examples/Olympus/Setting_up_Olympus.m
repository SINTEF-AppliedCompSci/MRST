mrstModule add deckformat ad-core

%%  Reading the deck for realization 1 of Olympus

fn = 'olympus1/OLYMPUS_1/OLYMPUS_1.DATA'
deck = readEclipseDeck(fn);
deck = convertDeckUnits(deck);
 
%%  Creating the MRST structures for simulation
[state0, model, schedule, nonlinear] = initEclipseProblemAD(deck, 'useMex', true);
problem = packSimulationProblem(state0, model, schedule, 'Olympus_model', 'NonLinearSolver', nonlinear);

%% Run simulation
[ok, status] = simulatePackedProblem(problem);

%% Getting simulation Results
 [wellSols, states, reports] = getPackedSimulatorOutput(problem);
    
 %% To run flow diagnostics vizualizer
 %d = PostProcessDiagnosticsMRST(problem);