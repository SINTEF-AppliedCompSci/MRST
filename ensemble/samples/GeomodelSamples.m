classdef GeomodelSamples < RockSamples
    
    methods
        %-----------------------------------------------------------------%
        function problem = setSample(samples, sampleData, problem)
            model  = problem.SimulatorSetup.model;
            model.rock = sampleData.rock;
            model.G = sampleData.G;
            model = model.setupOperators();
            problem.SimulatorSetup.model = model;
            W = problem.SimulatorSetup.schedule.control(1).W;
            for i = 1:max(model.G.cells.tag)
                W(i).cells = find(model.G.cells.tag == i);
            end
            problem.SimulatorSetup.schedule.control(1).W = W;
            if samples.updateWells
                % Update well indices
                schedule = problem.SimulatorSetup.schedule;
                schedule = samples.updateWellIndices(model, schedule);
                problem.SimulatorSetup.schedule = schedule;
            end
            state0 = problem.SimulatorSetup.state0;
            state0 = initResSol(model.G, state0.pressure(1), state0.s(1,:));
            problem.SimulatorSetup.state0 = state0;
        end
    end
    
end