classdef AutoDiffBackend
    % Automatic differentiation backend class
    %
    % SYNOPSIS:
    %   backend = AutoDiffBackend()
    %
    % DESCRIPTION:
    %   Base class for automatic differentiation backends. A backend is an
    %   implementation of automatic differentation and must support the
    %   initialization of one (or more) variables as AD objects.
    %
    % RETURNS:
    %  Backend - Initialized class instance
    %
    % SEE ALSO:
    %   `DiagonalAutoDiffBackend`, `SparseAutoDiffBackend`

    methods
        function backend = AutoDiffBackend()
            % Class constructor.
        end
        
        function model = updateDiscreteOperators(backend, model)
            % Modify model/operators for use with backend
            %
            % SYNOPSIS:
            %   model = backend.updateDiscreteOperators(model)
            %
            % DESCRIPTION:
            %   Detailed description of function
            %
            % PARAMETERS:
            %   backend - Class instance
            %
            % RETURNS:
            %   model - Model with updated operators
            %
            % NOTE:
            %   Base class implements no modifications. This function is
            %   automatically called as a part of `validateModel`.

        end
        
        function out = getBackendDescription(backend)
            % Get a text string describing the backend
            out = 'Standard ADI (sparse jacobian)';
        end
        
        function varargout = initVariablesAD(backend, varargin)
            % Initialize variables as AD
            n         = nargout;
            varargout = cell([1, n]);
            [varargout{:}] = initVariablesADI(varargin{1:n});
        end
        
        function v = convertToAD(backend, v, sample)
            % Given a `sample` AD object, convert `v` to AD.
            v = double2ADI(v, sample);
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
