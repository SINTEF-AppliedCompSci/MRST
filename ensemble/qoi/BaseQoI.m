classdef BaseQoI
    % Template class for extracting a quantity of interest from a simulated
    % problem.
    
    properties
        ResultHandler % Handler for writing/reading QoIs to/from file
    end
    
    methods
        
        %-----------------------------------------------------------------%
        function qoi = BaseQoI()
            % Constructor is intentionally empty
        end
        
        %-----------------------------------------------------------------%
        function qoi = validateQoI(qoi, problem)
            % Validate the quantity of interest. BaseQoI sets up an
            % appropriate ResultHandler, so that all subclass
            % implementations of this function should start with
            % qoi = validateQoI@BaseQoI(qoi, problem);
            if isempty(qoi.ResultHandler)
                % Set up ResultHandler. By default, data is stored in the
                % the dataDirectory of the problem output handler
                dataDir = problem.OutputHandlers.states.dataDirectory;
                qoi.ResultHandler = ResultHandler('dataDirectory', dataDir, ...
                                                  'dataFolder'   , ''     , ...
                                                  'dataPrefix'   , 'qoi' ); %#ok
            end
            % Check that output is stored with the correct name
            assert(strcmp(qoi.ResultHandler.dataPrefix, 'qoi'), ...
                   'ResultHandler data prefix must be ''qoi''.');
        end
        
        %-----------------------------------------------------------------%
        function u = getQoI(qoi, problem)
            % Get quantity of interest for a given problem
            % TODO: Check that problem has been simulated successfully, and
            % issue a warning if it is not
            seed = str2double(problem.OutputHandlers.states.dataFolder);
            if qoi.isComputed(seed)
                % QoI already computed - read from file
                u = qoi.ResultHandler{seed};
            else
                % Compute QoI and store to file
                u  = qoi.computeQoI(problem);
                us = u; % Handle special case when u is a cell array
                if iscell(us), us = {us}; end 
                qoi.ResultHandler{seed} = us;
            end 
        end
        
        %-----------------------------------------------------------------%
        function ok = isComputed(qoi, seed)
            % Check if qoi for a given seed it computed
            ids = qoi.ResultHandler.getValidIds();
            ok  = any(ids == seed);
        end
        
        %-----------------------------------------------------------------%
        function u = computeQoI(qoi, problem) %#ok
            % Compute quantity of interest for a given problem
            error('Template class not meant for direct use!');
        end
        
        %-----------------------------------------------------------------%
        function n = norm(qoi, u) %#ok
            % Compute norm of the quantity of interest
            n = abs(u);
        end
        
        %-----------------------------------------------------------------%
        function plotEnsemble(qoi, ensemble)
            % Create a meaningful plot of the ensemble based
            % on the relevant QoI.
            error('Template class not meant for direct use!');
        end
        
        %-----------------------------------------------------------------%
        function getEnsembleMean(qoi, ensemble)
            % Compute the mean quantity of interest for the given ensemble
            error('Template class not meant for direct use!');
        end
            
    end
end
    
