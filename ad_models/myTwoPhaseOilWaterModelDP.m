classdef myTwoPhaseOilWaterModelDP < TwoPhaseOilWaterModelDP
% Two phase oil/water system without dissolution
    properties

    end

    methods
        function model = myTwoPhaseOilWaterModelDP(G, rock, fluid, varargin)
            model = model@TwoPhaseOilWaterModelDP(G, rock, fluid, varargin{:});
            model = merge_options(model, varargin{:});
        end

        %--------------------------------------------------------------------%  
        function state = validateState(model, state)
            % Dummy validation to be called in simulateScheduleAD.
            return
        end             
        
    end
end

%{
Copyright 2009-2018 SINTEF Digital, Mathematics & Cybernetics.

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