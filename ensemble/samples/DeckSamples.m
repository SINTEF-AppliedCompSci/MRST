classdef DeckSamples < BaseSamples
    
    properties
        gridFromDeck = true % Construct MRST grid from deck
        initArgs     = {}   % Arguments passed directly to initEclipseProblemAD
    end
    
    methods
        function problem = setSample(samples, sampleData, problem)
            if samples.gridFromDeck
                % Construct MRST grid from deck
                G = [];
            else
                % Use grid from problem
                G = problem.simulatorSetup.model.G;
            end
            % Get initial state, model and schedule
            [state0, model, schedule] = initEclipseProblemAD(sampleData, ...
                                                 'G', G, samples.initArgs{:});
            % update problem
            problem.SimulatorSetup.state0   = state0;
            problem.SimulatorSetup.model    = model;
            problem.SimulatorSetup.schedule = schedule;
        end
    end
    
end