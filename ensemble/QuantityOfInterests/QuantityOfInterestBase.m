classdef QuantityOfInterestBase
    % Template class for extracting a quantity of interest from a simulated
    % problem.
    
    properties
        % No properties
    end
    
    methods
        
        function quantityOfInterestBase = QuantityOfInterestBase()
            % Constructor is intentionally empty
        end
        
        function ok = save(quantityOfInterestBase, path, qoi)
            error('Implementation missing');
        end
        
        function qoi = getQoI(quantityOfInterestBase, problem)
            % getQoI reads the result files of the given problem and
            % extracts the quantity of interest, which is returned.
            error('Template class not meant for direct use!');
        end
        
        function n = norm(quantityOfInterestBase, qoi)
            n = abs(qoi);
        end
    end
end
    
