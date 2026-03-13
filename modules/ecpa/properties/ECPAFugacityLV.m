classdef ECPAFugacityLV < StateFunction
    % Fugacity for liquid-vapor system
    properties
        useCompactEvaluation = true;
    end
    
    methods
        function gp = ECPAFugacityLV(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'PhaseMixingCoefficients', 'ComponentPhaseMoleFractions', 'PhaseCompressibilityFactors'});
            gp = gp.dependsOn({'pressure', 'T'}, 'state');
            gp.label = 'f_\alpha';
        end

        function f = evaluateOnDomain(prop, model, state)
            eos = model.EOSModel;
            [p, T] = model.getProps(state, 'pressure', 'T');
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
                z = Z{phix};
                if i == 2 && prop.useCompactEvaluation && ~all(twoPhase)
                    [~, ~, twoPhase] = model.getFlag(state);
                    fi = f(:, 1);
                    if any(twoPhase)
                        if iscell(xy)
                            xy = cellfun(@(x) x(twoPhase), xy, 'UniformOutput', false);
                        else
                            xy = xy(twoPhase, :);
                        end
                        if iscell(m.a)
                            Ai = cellfun(@(x) x(twoPhase), m.Ai, 'UniformOutput', false);
                        else
                            Ai = Ai(twoPhase, :);
                        end
                        if iscell(m.acti)
                            Si = cellfun(@(x) x(twoPhase), m.Si, 'UniformOutput', false);
                        else
                            Si = Si(twoPhase, :);
                        end
                        if iscell(m.Tr)
                            Tr = cellfun(@(x) x(twoPhase), m.Tr, 'UniformOutput', false);
                        else
                            Tr = Tr(twoPhase, :);
                        end
                        if iscell(z.XA)
                            XA = cellfun(@(x) x(twoPhase), z.XA, 'UniformOutput', false);
                        else
                            XA = z.XA(twoPhase, :);
                        end
                        f_2ph = model.EOSModel.computeFugacity(T(twoPhase),p(twoPhase),xy,m.A(twoPhase), m.B(twoPhase),Ai, z.v(twoPhase),Si, XA, z.XC(twoPhase), Tr);
                        for j = 1:numel(fi)
                            fi{j}(twoPhase) = f_2ph{j};
                        end
                    end
                    f(:, i) = fi;
                else
                    f(:, i) = model.EOSModel.computeFugacity(T, p, xy, m.A, m.B, m.Ai, z.v, m.Si, z.XA, z.XC, m.Tr);
                end
            end
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
