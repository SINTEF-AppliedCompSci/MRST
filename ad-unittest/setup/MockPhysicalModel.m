classdef MockPhysicalModel < PhysicalModel
    
    methods
        function model = MockPhysicalModel()
            model = model@PhysicalModel([]);
        end
        
        function [fn, index] = getVariableField(model, name, varargin)
            % Get the index/name mapping for the model (such as where
            % pressure or water saturation is located in state)
            switch(lower(name))
                case 'x1'
                    fn = 'x';
                    index = 1;
                case 'x2'
                    fn = 'x';
                    index = 2;
                case 'x'
                    fn = 'x';
                    index = ':';
                case 'y'
                    fn = 'y';
                    index = 1;
                otherwise
                    % This will throw an error for us
                    [fn, index] = getVariableField@PhysicalModel(model, name, varargin{:});
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
