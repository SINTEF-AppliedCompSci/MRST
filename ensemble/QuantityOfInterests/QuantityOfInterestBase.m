classdef QuantityOfInterestBase
    % Template class for extracting a quantity of interest from a simulated
    % problem.
    
    properties
        % No properties
    end
    
    methods
        
        function quantityOfInterestBase = QuantityOfInterestBase()
            %error('Template class not meant for direct use!');
            % Some silly comment just to get another checksum and more
        end
        
        function ok = save(quantityOfInterestBase, path, qoi)
            error('Implementation missing');
        end
        
        function qoi = getQOI(quantityOfInterestBase, problem)
            % getQOI reads the result files of the given problem and
            % extracts the quantity of interest, which is returned.
            error('Template class not meant for direct use!');
        end
        
        function n = norm(quantityOfInterestBase, qoi)
            n = abs(qoi);
        end
    end
end
    
