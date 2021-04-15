classdef SparseAutoDiffBackend < AutoDiffBackend
    % Automatic differentiation backend class (sparse representation)
    %
    % SYNOPSIS:
    %   backend = SparseAutoDiffBackend()
    %
    % DESCRIPTION:
    %   This version of the AD backend uses different types of sparse
    %   blocks to represent derivatives.
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
            out = 'Sparse';
            if ~backend.useBlocks
                out = [out, ' (single matrix)'];
            end
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
                v = double2GenericAD(v, sample);
            end
        end
    end
end

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
