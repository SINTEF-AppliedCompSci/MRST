classdef SparseAutoDiffBackend < AutoDiffBackend
    % Automatic differentiation backend class (sparse representation)
    %
    % SYNOPSIS:
    %   backend = SparseAutoDiffBackend()
    %
    % DESCRIPTION:
    %    This version of the AD backend uses different types of sparse
    %    blocks to represent derivatives.
    %
    % RETURNS:
    %   Backend - Initialized class instance
    %
    % SEE ALSO:
    %   `AutoDiffBackend`, `DiagonalAutoDiffBackend`

    properties
        useBlocks % Organize Jacobian as a set of blocks, instead of one large AD matrix
    end
    
    methods
        function backend = SparseAutoDiffBackend(varargin)
            backend = backend@AutoDiffBackend();
            backend.useBlocks = true;
            backend = merge_options(backend, varargin{:});
        end
        
        function model = updateDiscreteOperators(backend, model)
            % Do nothing
        end
        
        function out = getBackendDescription(backend)
            out = 'Standard ADI (sparse jacobian)';
        end
        
        function varargout = initVariablesAD(backend, varargin)
            n         = nargout;
            varargout = cell([1, n]);
            if backend.useBlocks
                [varargout{:}] = initVariablesADI(varargin{1:n});
            else
                [varargout{:}] = initVariablesAD_oneBlock(varargin{1:n});
            end
        end
        
        function v = convertToAD(backend, v, sample)
            if backend.useBlocks
                v = double2ADI(v, sample);
            else
                v = double2NewAD(v, sample);
            end
        end
        
    end
end