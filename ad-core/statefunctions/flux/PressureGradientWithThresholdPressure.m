classdef PressureGradientWithThresholdPressure < PressureGradient
    properties
        pressureThreshold
    end
    
    methods
        function gp = PressureGradientWithThresholdPressure(model, varargin)
            gp@PressureGradient(model, varargin{:});
            equil = model.inputdata.REGIONS.EQLNUM(model.G.cells.indexMap);
            thpres = model.inputdata.SOLUTION.THPRES;
            N = model.operators.N;
            nf = size(N, 1);
            threshold = zeros(nf, 2);
            
            left_equil = equil(N(:, 1));
            right_equil = equil(N(:, 2));
            for i = 1:size(thpres, 1)
                left = left_equil == thpres(i, 1);
                right = right_equil == thpres(i, 2);
                
                threshold(left & right, 1) = thpres(i, 3);
                % Other direction
                left = left_equil == thpres(i, 2);
                right = right_equil == thpres(i, 1);
                threshold(left & right, 2) = thpres(i, 3);
            end
            gp.pressureThreshold = threshold;
        end
        function dp = evaluateOnDomain(prop, model, state)
            dp = evaluateOnDomain@PressureGradient(prop, model, state);
            p = model.getProps(state, 'pressure');
            pv = value(p);
            left = model.operators.N(:, 1);
            right = model.operators.N(:, 2);
            pL = pv(left);
            pR = pv(right);
            
            dpv = pL - pR;
            % Flow from left to right
            posFlow = pL > pR;
            thres = prop.pressureThreshold;
            act = dpv > thres(:, 1) | dpv < thres(:, 2);
            delta = thres(:, 2);
            delta(posFlow) = thres(posFlow, 1);
            for i = 1:numel(dp)
                dp{i} = act.*(dp{i} - delta);
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
