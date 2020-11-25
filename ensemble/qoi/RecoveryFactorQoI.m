classdef RecoveryFactorQoI < BaseQoI
    
    properties
        phase = 'oil';
    end
    
    methods
        
        %-----------------------------------------------------------------%
        function qoi = validateQoI(qoi, problem, varargin)
            % Let base class do its thing
            qoi   = validateQoI@BaseQoI(qoi, problem, varargin{:});
            % Check that phase exists
            model = problem.SimulatorSetup.model;
            [~, phases] = model.getPhaseNames();
            assert(any(strcmpi(qoi.phase, phases)), ...
                   'Did not find %s in model phases', qoi.phase);
        end
        
        %-----------------------------------------------------------------%
        function u = computeQoI(qoi, problem)
            % Get state output handler
            [~, states] = getPackedSimulatorOutput(problem, 'readFromDisk', false);
            % Get model
            model = problem.SimulatorSetup.model.validateModel();
            % Compute initial volume
            [si, pvi] = model.getProps(states{1}, qoi.phase, 'PoreVolume');
            vi = sum(si.*pvi);
            % Compute final volume
            [sf, pvf] = model.getProps(states{numelData(states)}, qoi.phase, 'PoreVolume');
            vf = sum(sf.*pvf);
            % Compute recovery factor
            u = {(vi-vf)./vi};
        end
        
    end

end