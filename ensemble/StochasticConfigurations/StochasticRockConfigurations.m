classdef StochasticRockConfigurations < StochasticConfigurations
    
    properties
        %
        % 
    end
    
    methods
        
        function stochRockConfig = StochasticRockConfigurations(data)
            stochRockConfig = stochRockConfig@StochasticConfigurations();
            
            stochRockConfig.data = data;
            if ~isempty(data)
                stochRockConfig.num = numel(data);
            end
        end
            
        
        function specificProblem = getProblem(stochRockConfig, baseProblem, seed)
            assert(seed < stochRockConfig.num, ...
                'seed is larger than number of ensemble members');
            
            if isinf(stochRockConfig.num)
                rockConfig = stochRockConfig.sampleConfig(seed);
            else
                rockConfig = stochRockConfig.data{seed};
            end
            
            newRock = makeRock(baseProblem.SimulatorSetup.model.G, ...
                rockConfig.perm(:), rockConfig.poro(:));
            
            specificProblem = baseProblem;
            specificProblem.SimulatorSetup.model.rock = newRock;
            
            specificProblem.SimulatorSetup.model.operators ...
                = setupOperatorsTPFA(problem.SimulatorSetup.model.G, ...
                    problem.SimulatorSetup.model.rock);
            
            specificProblem.SimulatorSetup.model = problem.SimulatorSetup.model.removeStateFunctionGroupings();
            
        end        
        
        
    end
end