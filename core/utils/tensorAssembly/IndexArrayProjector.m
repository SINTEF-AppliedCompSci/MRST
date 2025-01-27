classdef IndexArrayProjector 

    properties
        
        tbl        % IndexArray we project
        projfds    % fields that we use for the projection
        replacefds % Possibility to change field names in table
                   % (before setting projection tbl)

        %% dispactch mappings updated during setup
        inds
        
    end

    methods

        function proj = IndexArrayProjector(varargin)

            proj.tbl        = [];
            proj.replacefds = [];
            proj.projfds    = [];

            proj = merge_options(proj, varargin{:}); 

            if isempty(proj.projfds)
                proj.projfds = {};
            end            

        end

        function [tbl, proj] = eval(proj)

            projfds = proj.projfds;
            tbl     = proj.tbl;

            if ~isempty(proj.replacefds)
                 tbl = replacefield(tbl, proj.replacefds);
            end

            [tbl, indstruct] = projIndexArray(tbl, projfds);
            
            proj.tbl  = tbl;
            proj.inds = indstruct.inds;
            
        end

    end
    
end

%{
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

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


