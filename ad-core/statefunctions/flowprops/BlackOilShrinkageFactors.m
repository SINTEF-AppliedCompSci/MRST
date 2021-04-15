classdef BlackOilShrinkageFactors < ShrinkageFactors
    % Shrinkage factors for black-oil
    properties
        useSaturatedFlag = true;
        disgas = false;
        vapoil = false;
    end
    
    methods
        function b = BlackOilShrinkageFactors(model, varargin)
            b@ShrinkageFactors(model, varargin{:});
            assert(isa(model, 'ThreePhaseBlackOilModel'), ...
                ['Model must be derived from the black-oil model. Did you want', ...
                ' the regular ShrinkageFactors class instead?'])
            b.disgas = model.disgas;
            b.vapoil = model.vapoil;
            if b.disgas
                b = b.dependsOn({'rs'}, 'state');
            end
            if b.vapoil
                b = b.dependsOn({'rv'}, 'state');
            end
            if b.useSaturatedFlag
                b = b.dependsOn({'s'}, 'state');
            end
        end

        function b = evaluatePhaseShrinkageFactor(prop, model, state, name, p)
            if prop.disgas && strcmp(name, 'O')
                % Oileic phase with dissolved gas component
                rs = model.getProp(state, 'rs');
                if prop.useSaturatedFlag
                    sG = model.getProp(state, 'sg');
                    flag = sG > 0;
                else
                    flag = false(numelValue(rs), 1);
                end
                b = prop.evaluateFluid(model, 'bO', p, rs, flag);
            elseif prop.vapoil && strcmp(name, 'G')
                % Gaseous phase with vaporized oil component
                rv = model.getProp(state, 'rv');
                if prop.useSaturatedFlag
                    sO = model.getProp(state, 'so');
                    flag = sO > 0;
                else
                    flag = false(numelValue(rv), 1);
                end
                b = prop.evaluateFluid(model, 'bG', p, rv, flag);
            else
                % Can use base class directly
                b = evaluatePhaseShrinkageFactor@ShrinkageFactors(prop, model, state, name, p);
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
