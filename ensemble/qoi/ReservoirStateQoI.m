classdef ReservoirStateQoI < BaseQoI
    
    properties
        name = 'sW'
        cells = ':'
        time  = inf
        
        dt = nan
        dx
    end
    
    methods
        %-----------------------------------------------------------------%
        function qoi = ReservoirStateQoI(varargin)
            qoi = qoi@BaseQoI();
            qoi = merge_options(qoi, varargin{:});
        end
        
        %-----------------------------------------------------------------%
        function qoi = validateQoI(qoi, problem)
            qoi = validateQoI@BaseQoI(qoi, problem);
            qoi.dx = problem.SimulatorSetup.model.G.cells.volumes;
            if strcmpi(qoi.time, ':')
                qoi.time = cumsum(problem.SimulatorSetup.schedule.step.val);
            end
            if numel(qoi.time) == 1
                if isinf(qoi.time)
                    qoi.time = sum(problem.SimulatorSetup.schedule.step.val);
                end
            else
                qoi.dt = diff([0,qoi.time]);
            end
        end
        
        %-----------------------------------------------------------------%
        function u = computeQoI(qoi, problem)
            nt = numel(qoi.time);
            u  = cell(nt,1);
            states   = problem.OutputHandlers.states;
            model    = problem.SimulatorSetup.model.setupStateFunctionGroupings();
            schedule = problem.SimulatorSetup.schedule;
            for i = 1:nt
                n     = qoi.findTimestep(schedule, qoi.time(i));
                state = states{n};
                u{i}  = qoi.getStateValue(model, state);
            end
        end
        
        %-----------------------------------------------------------------%
        function n = findTimestep(qoi, schedule, time)
            [~, n] = min(abs(cumsum(schedule.step.val) - time));
        end
        
        %-----------------------------------------------------------------%
        function u = getStateValue(qoi, model, state)
            u = model.getProp(state, qoi.name);
        end
        
        %-----------------------------------------------------------------%
        function n = norm(qoi, u)
            n = cellfun(@(u) sum(abs(u).*qoi.dx), u);
            if ~isnan(qoi.dt)
                n = n.*qoi.dt;
            end
            n = sum(n);
        end
        
    end
    
end