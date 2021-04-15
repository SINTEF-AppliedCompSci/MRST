classdef ComponentPhaseMoleFractionsLV < StateFunction
    % Mole fractions for each phase
    properties
    end
    
    methods
        function gp = ComponentPhaseMoleFractionsLV(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'x', 'y'}, 'state');
            gp.label = 'x_{i,\alpha}';
        end

        function v = evaluateOnDomain(prop, model, state)
            [x, y] = model.getProps(state, 'x', 'y');
            if ~iscell(x)
                x = expandMatrixToCell(x);
            end
            if ~iscell(y)
                y = expandMatrixToCell(y);
            end
            % Ncomp by nph matrix. The eos components get the two-phase
            % distribution and the immiscible components get the identity
            nph = model.getNumberOfPhases();
            ncomp = model.getNumberOfComponents();
            nc = model.G.cells.num;
            v = cell(ncomp, nph);
            isEoS = model.getEoSComponentMask();
            [v{isEoS, model.getLiquidIndex()}] = x{:};
            [v{isEoS, model.getVaporIndex()}] = y{:};
            if ~all(isEoS)
                if isa(model, 'GenericReservoirModel')
                    for cNo = find(~isEoS)
                        ix = model.Components{cNo}.phaseIndex;
                        v{cNo, ix} = ones(nc, 1);
                    end
                else
                    % Just water - first component
                    v{1, 1} = ones(nc, 1);
                end
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
