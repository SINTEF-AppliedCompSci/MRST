classdef ComponentPhaseMoleFractionsLV < StateFunction
    properties
    end
    
    methods
        function gp = ComponentPhaseMoleFractionsLV(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'x', 'y'}, 'state');
        end

        function v = evaluateOnDomain(prop, model, state)
            [x, y] = model.getProps(state, 'x', 'y');
            if ~iscell(x)
                x = expandMatrixToCell(x);
            end
            if ~iscell(y)
                y = expandMatrixToCell(y);
            end
            if model.water
                ncomp = numel(x) + 1;
                nc = model.G.cells.num;
                u = ones(nc, 1);
                x = [{[]}, x];
                y = [{[]}, y];
                w = cell(1, ncomp);
                w{1} = u;
                v = [w', x', y'];
            else
                v = [x', y'];
            end
        end
    end
end

%{
Copyright 2009-2019 SINTEF Digital, Mathematics & Cybernetics.

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
