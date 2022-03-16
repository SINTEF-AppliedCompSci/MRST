classdef MCLevelQoI < BaseQoI
    
    properties
        levelQoIs
    end
    
    methods
        function qoi = MCLevelQoI(levelQoIs)
            cls = cellfun(@class, levelQoIs, 'UniformOutput', false);
            assert(all(strcmpi(cls{1}, cls)));
            assert(numel(levelQoIs) <= 2);
            qoi = qoi@BaseQoI();
            qoi.levelQoIs = levelQoIs;
            qoi.names = qoi.levelQoIs{1}.names;
        end
        
        function du = getQoI(qoi, seed)
            du = qoi.levelQoIs{end}.ResultHandler{seed};
            if numel(qoi.levelQoIs) == 2
                u0 = qoi.levelQoIs{1}.ResultHandler{seed};
                for fn = qoi.names
                    du.(fn{1}) = du.(fn{1}) - u0.(fn{1});
                end
                du.cost = du.cost + u0.cost;
            end
            qoi.ResultHandler{seed} = du;
        end
        
    end
end

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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
