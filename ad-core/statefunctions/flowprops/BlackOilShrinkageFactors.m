classdef BlackOilShrinkageFactors < StateFunction
    % Shrinkage factors for black-oil
    properties
        useSaturatedFlag = true;
        disgas = false;
        vapoil = false;
    end
    
    methods
        function gp = BlackOilShrinkageFactors(model, varargin)
            gp@StateFunction(model, varargin{:});
            if isprop(model, 'disgas')
                gp.disgas = model.disgas;
                if gp.disgas
                    gp = gp.dependsOn({'rs'}, 'state');
                end
            end
            if isprop(model, 'vapoil')
                gp.vapoil = model.vapoil;
                if gp.vapoil
                    gp = gp.dependsOn({'rv'}, 'state');
                end
            end
            gp = gp.dependsOn({'PhasePressures'});
            gp.label = 'b_\alpha';
        end

        function b = evaluateOnDomain(prop, model, state)
            [act, phInd] = model.getActivePhases();
            nph = sum(act);
            b = cell(1, nph);
            
            f = model.fluid;
            p_phase = prop.getEvaluatedDependencies(state, 'PhasePressures');
            if model.water
                wix = phInd == 1;
                pw = p_phase{wix};
                bW = prop.evaluateFunctionOnDomainWithArguments(f.bW, pw);
                b{wix} = bW;
            end
            
            if model.oil
                oix = phInd == 2;
                po = p_phase{oix};
                if prop.disgas
                    rs = model.getProp(state, 'rs');
                    if prop.useSaturatedFlag
                        sG = model.getProp(state, 'sg');
                        flag = sG > 0;
                    else
                        flag = false(numelValue(po), 1);
                    end
                    bO = prop.evaluateFunctionOnDomainWithArguments(f.bO, po, rs, flag);
                else
                    bO = prop.evaluateFunctionOnDomainWithArguments(f.bO, po);
                end
                b{oix} = bO;
            end
            
            if model.gas
                gix = phInd == 3;
                pg = p_phase{gix};
                if prop.vapoil
                    rv = model.getProp(state, 'rv');
                    if prop.useSaturatedFlag
                        sO = model.getProp(state, 'so');
                        flag = sO > 0;
                    else
                        flag = false(numelValue(pg), 1);
                    end
                    bG = prop.evaluateFunctionOnDomainWithArguments(f.bG, pg, rv, flag);
                else
                    bG = prop.evaluateFunctionOnDomainWithArguments(f.bG, pg);
                end
                b{gix} = bG;
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
