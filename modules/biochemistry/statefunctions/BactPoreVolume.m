classdef BactPoreVolume < StateFunction
    % Static pore volume accessor for bacterial modeling
    properties
        % No additional properties needed
    end
    
    methods
        function gp = BactPoreVolume(model, varargin)
            % Constructor with pore volume validation
            gp@StateFunction(model, varargin{:});
            
            % Validate pore volume field exists and is properly sized
            assert(isfield(model.operators, 'pv'), ...
                'Pore-volume (pv) must be present as field in operators struct');
            assert(numelValue(model.operators.pv) == model.G.cells.num,...
                'Pore-volumes must be defined in each cell.');
            assert(all(model.operators.pv > 0), ...
                'Pore-volumes must be non-negative.');
            
            gp.label = 'Phi';
        end
        
        function pv = evaluateOnDomain(prop, model, state) %#ok
            % Retrieve static pore volume values
            %
            % PARAMETERS:
            %   prop  - Property function instance
            %   model - Reservoir model instance
            %   state - State struct (unused for static PV)
            %
            % RETURNS:
            %   pv    - Pore volume values for all cells
            
            % Direct access to pre-computed pore volumes
            pv = model.operators.pv;
        end
    end
end

%{
Copyright 2009-2023 SINTEF Digital, Mathematics & Cybernetics.

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