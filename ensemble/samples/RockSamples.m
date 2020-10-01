classdef RockSamples < BaseSamples

    methods
        function problem = setSample(samples, sampleData, problem)
            warning('Assuming fully implicit model');
            % TODO: Add support for derived models (SI, Wrappers, ...)
            % Get model
            model = problem.SimulatorSetup.model;
            % Set rock properties from sample data
            model.rock.perm = sampleData.perm;
            model.rock.poro = sampleData.poro;
            % Set up operators
            model = model.setupOperators();
            % Replace problem model
            problem.SimulatorSetup.model = model;
        end
    end
    
end