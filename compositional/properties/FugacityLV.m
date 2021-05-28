classdef FugacityLV < StateFunction
    % Fugacity for liquid-vapor system
    properties
        useCompactEvaluation = true;
    end
    
    methods
        function gp = FugacityLV(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'PhaseMixingCoefficients', 'ComponentPhaseMoleFractions', 'PhaseCompressibilityFactors'});
            gp = gp.dependsOn({'pressure'}, 'state');
            gp.label = 'f_\alpha';
        end

        function f = evaluateOnDomain(prop, model, state)
            eos = model.EOSModel;
            p = model.getProps(state, 'pressure');
            [mix, mf, Z] = prop.getEvaluatedDependencies(state, 'PhaseMixingCoefficients', 'ComponentPhaseMoleFractions', 'PhaseCompressibilityFactors');
            
            ncomp = eos.getNumberOfComponents();
            f = cell(ncomp, 2);
            twoPhase = model.getTwoPhaseFlag(state);
            isEoS = model.getEoSComponentMask();

            for i = 1:2
                if i == 1
                    phix = model.getLiquidIndex();
                else
                    phix = model.getVaporIndex();
                end
                xy = mf(isEoS, phix)';
                m = mix{phix};
                if i == 2 && prop.useCompactEvaluation && ~all(twoPhase)
                    [~, ~, twoPhase] = model.getFlag(state);
                    fi = f(:, 1);
                    if any(twoPhase)
                        xy = cellfun(@(x) x(twoPhase), xy, 'UniformOutput', false);
                        Si = cellfun(@(x) x(twoPhase), m.Si, 'UniformOutput', false);
                        Bi = cellfun(@(x) x(twoPhase), m.Bi, 'UniformOutput', false);
                        f_2ph = model.EOSModel.computeFugacity(p(twoPhase), xy, Z{phix}(twoPhase), m.A(twoPhase), m.B(twoPhase), Si, Bi);
                        for j = 1:numel(fi)
                            fi{j}(twoPhase) = f_2ph{j};
                        end
                    end
                    f(:, i) = fi;
                else
                    f(:, i) = model.EOSModel.computeFugacity(p, xy, Z{phix}, m.A, m.B, m.Si, m.Bi);
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
