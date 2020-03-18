classdef BlackOilViscosity < StateFunction
    % Black-oil style viscosity functions that account for rs and Rv
    properties
        useSaturatedFlag = true;
        disgas = false;
        vapoil = false;
    end
    
    methods
        function gp = BlackOilViscosity(model, varargin)
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
            gp.label = '\mu_\alpha';
        end
        
        function mu = evaluateOnDomain(prop, model, state)
            [act, phInd] = model.getActivePhases();
            nph = sum(act);
            mu = cell(1, nph);
            
            f = model.fluid;
            p_phase = prop.getEvaluatedDependencies(state, 'PhasePressures');
            [sample, isAD] = getSampleAD(p_phase{:});
            nc = numelValue(sample);

            if model.water
                wix = phInd == 1;
                pw = p_phase{wix};
                mu{wix} = prop.evaluateFluid(model, 'muW', pw);
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
                        flag = false(nc, 1);
                    end
                    extra = {rs, flag};
                else
                    extra = {};
                end
                mu{oix} = prop.evaluateFluid(model, 'muO', po, extra{:});
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
                        flag = false(nc, 1);
                    end
                    extra = {rv, flag};
                else
                    extra = {};
                end
                mu{gix} = prop.evaluateFluid(model, 'muG', pg, extra{:});
            end
            if isAD
                for i = 1:numel(mu)
                    if ~isa(mu{i}, 'ADI')
                        mu{i} = model.AutoDiffBackend.convertToAD(mu{i}, sample);
                    end
                end
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
