classdef AutoDiffBackend
    properties
        
    end
    
    methods
        function backend = AutoDiffBackend()
            
        end
        
        function ops = updateDiscreteOperators(backend, ops)
            % Do nothing
        end
        
        function out = getBackendDescription(backend)
            out = 'Standard ADI (sparse jacobian)';
        end
        
        function varargout = initVariablesAD(backend, varargin)
            n         = nargout;
            varargout = cell([1, n]);
            [varargout{:}] = initVariablesADI(varargin{1:n});
        end
        
        function v = convertToAD(backend, v, sample)
            v = double2ADI(v, sample);
        end
    end
end