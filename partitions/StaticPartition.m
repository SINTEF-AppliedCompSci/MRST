classdef StaticPartition < Partition
    % Static partition based on precomputed partition vector
    
    properties
    end
    
    methods
        %-----------------------------------------------------------------%
        function sp = StaticPartition(partition)
            sp = sp@Partition('value', partition);
            sp.numBlocks = max(partition);
        end
        
        %-----------------------------------------------------------------%
        function ok = doCompute(p, model, state, state0, dt, drivingForces) %#ok
            ok = false;
        end
        
    end
    
end

%{
Copyright 2009-2020 SINTEF Digital, Mathematics & Cybernetics.

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
