classdef PhasePotentialDifferenceThresholded < PhasePotentialDifference
    properties
        pressureThreshold
    end
    
    methods
        function gp = PhasePotentialDifferenceThresholded(model, varargin)
            gp@PhasePotentialDifference(model, varargin{:});
            deck = model.inputdata;
            equil = deck.REGIONS.EQLNUM(model.G.cells.indexMap);
            thpres = deck.SOLUTION.THPRES;
            irreversible = isfield(deck.RUNSPEC, 'EQLOPTS') && any(strcmp(deck.RUNSPEC.EQLOPTS, 'IRREVERS'));
            if ~irreversible
                % Repeat twice to get both directions
                thpres = [thpres; thpres(:, [2, 1, 3])];
            end

            N = model.operators.N;
            nf = size(N, 1);
            threshold = zeros(nf, 2);
            
            left_equil = equil(N(:, 1));
            right_equil = equil(N(:, 2));
            for i = 1:size(thpres, 1)
                left = left_equil == thpres(i, 1);
                right = right_equil == thpres(i, 2);
                
                threshold(left & right, 1) = -thpres(i, 3);
                % Other direction
                left = left_equil == thpres(i, 2);
                right = right_equil == thpres(i, 1);
                threshold(left & right, 2) = thpres(i, 3);
            end
            if any(~isfinite(threshold(:)))
                warning('Non-defaulted threshold pressures detected. Please call setThresholdPressuresFromState');
            end
            gp.pressureThreshold = threshold;
        end

        function gp = setThresholdPressuresFromState(gp, model, state)
            % Defaulted values must be taken from initial conditions
            pt = gp.pressureThreshold;
            % Temp change to avoid recursion
            model.FlowDiscretization.PhasePotentialDifference = PhasePotentialDifference(model);
            state = model.validateState(state);
            state = model.getStateAD(state, true);
            
            pot = value(model.getProp(state, 'PhasePotentialDifference'));

            nf = size(pot, 1);
            s = value(state.s);
            T = true(nf, 1);
            % Negative potential - true upwind flag
            s_up = model.operators.faceUpstr(T, s);
            s_down = model.operators.faceUpstr(~T, s);

            % Ignore potential if phase is not present
            [pot_pos, pot_neg] = deal(pot);
            pot_pos(s_down == 0) = nan;
            pot_neg(s_up == 0) = nan;
            
            posflow = min(pot_neg, [], 2);
            negflow = max(pot_pos, [], 2);
            
            assert(all(isfinite(posflow)))
            assert(all(isfinite(negflow)))
            
            replace_pos = isnan(pt(:, 1)) & isfinite(posflow);
            replace_neg = isnan(pt(:, 2)) & isfinite(negflow);
            
            pt(replace_pos, 1) = posflow(replace_pos);
            pt(replace_neg, 2) = negflow(replace_neg);
            
            pt(isnan(pt)) = 0;
            gp.pressureThreshold = pt;
        end
        
        function pot = evaluateOnDomain(prop, model, state)
            pot = evaluateOnDomain@PhasePotentialDifference(prop, model, state);
            dpv = value(pot);
            thres = prop.pressureThreshold;
            
            for i = 1:numel(pot)
                dpi = dpv(:, i);
                posFlow = dpi < 0;
                
                delta = thres(:, 2);
                delta(posFlow) = thres(posFlow, 1);
                
                dp = pot{i};
                dp = dp - delta;
                ok = sign(dpi) == sign(value(dp));
                dp = dp.*ok;
                pot{i} = dp;
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
