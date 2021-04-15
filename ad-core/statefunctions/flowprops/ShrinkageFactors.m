classdef ShrinkageFactors < StateFunction
    % Shrinkage factors that depend only on pressure
    properties
        usePhasePressures = true;
    end
    
    methods
        function b = ShrinkageFactors(model, varargin)
            b@StateFunction(model, varargin{:});
            if b.usePhasePressures
                b = b.dependsOn('PhasePressures');
            else
                b = b.dependsOn('pressure', 'state');
            end
            b.label = 'b_\alpha';
            b.outputRange = [0, inf];
        end
        
        function b = evaluateOnDomain(prop, model, state)
            names = model.getPhaseNames();
            nph = numel(names);
            b = cell(1, nph);
            if prop.usePhasePressures
                p_phase = prop.getEvaluatedDependencies(state, 'PhasePressures');
            else
                p_phase = cell(1, nph);
                p = model.getProp(state, 'pressure');
                [p_phase{:}] = deal(p);
            end
            [sample, isAD] = getSampleAD(p_phase{:});
            for ph = 1:nph
                b{ph} = prop.evaluatePhaseShrinkageFactor(model, state, names(ph), p_phase{ph});
            end
            if isAD
                for i = 1:numel(b)
                    if ~isa(b{i}, 'ADI')
                        b{i} = model.AutoDiffBackend.convertToAD(b{i}, sample);
                    end
                end
            end
        end
        
        function mu = evaluatePhaseShrinkageFactor(prop, model, state, name, p) %#ok
            mu = prop.evaluateFluid(model, ['b', name], p);
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
