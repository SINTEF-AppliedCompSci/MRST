classdef BaseQoI
    % Template class for extracting a quantity of interest from a simulated
    % problem.
    
    properties
        % No properties
    end
    
    methods
        
        %-----------------------------------------------------------------%
        function qoi = BaseQoI()
            % Constructor is intentionally empty
        end
        
        %-----------------------------------------------------------------%
        function qoi = validateQoI(qoi, problem)
            % Function that potentially updates the object with properties
            % from the baseProblem.
            % This function is for instance called from in the 
            % MRSTEnsemble constructor.
            
            % Intentionally empty, as it in its simplest form just returns
            % itself.
        end
        
        %-----------------------------------------------------------------%
        function qoi = getQoI(qoi, problem)
            % getQoI reads the result files of the given problem and
            % extracts the quantity of interest, which is returned.
            error('Template class not meant for direct use!');
        end
        
        function n = norm(qoi, u)
            n = abs(u);
        end
    end
end
    
