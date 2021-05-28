classdef BlackOilPressureReductionFactors < PressureReductionFactors
    % Component weighting factors used to form a pressure equation
    properties
        disgas = false;
        vapoil = false;
        useUndersaturated = [];
    end
    
    methods
        function w_p = BlackOilPressureReductionFactors(model, varargin)
            w_p@PressureReductionFactors(model, varargin{:});
            assert(isa(model, 'ThreePhaseBlackOilModel'), ...
                ['Model must be derived from the black-oil model. ', ...
                'Did you want the regular PressureReductionFactors class instead?'])
            w_p = w_p.resetDependencies();
            w_p.vapoil = model.vapoil;
            w_p.disgas = model.disgas;
            if w_p.disgas
                w_p = w_p.dependsOn('rs', 'state');
                if ~w_p.vapoil && isempty(w_p.useUndersaturated)
                    w_p.useUndersaturated = true;
                end
            elseif isempty(w_p.useUndersaturated)
                w_p.useUndersaturated = false;
            end
            if w_p.vapoil
                w_p = w_p.dependsOn('rv', 'state');
            end
            w_p = w_p.dependsOn({'ShrinkageFactors', 'PoreVolume', 'Density'});
        end

        function w = evaluateOnDomain(prop, model, state)
            [b, pv, rho] = prop.getEvaluatedDependencies(state, 'ShrinkageFactors', 'PoreVolume', 'Density');
            rhoS = model.getSurfaceDensities();
            
            vap = prop.vapoil;
            dis = prop.disgas;
            if vap || dis
                oix = model.getPhaseIndex('O');
                gix = model.getPhaseIndex('G');
            end
            
            if vap
                rv = model.getProp(state, 'rv');
            else
                rv = 0;
            end
            
            if dis
                rs = model.getProp(state, 'rs');
                doUsat = prop.useUndersaturated;
                if doUsat
                    assert(~model.vapoil);
                    oix = model.getPhaseIndex('O');
                    bo = b{oix};
                    sg = model.getProp(state, 'sg');
                    p = model.getProp(state, 'pressure');
                    usat = value(sg) == 0;
                    rsAD = initVariablesAD_diagonal(value(rs));
                    tmp = model.PVTPropertyFunctions.ShrinkageFactors.evaluateFunctionOnDomainWithArguments(model.fluid.bO, value(p), rsAD, false(numelValue(rs), 1));
                    if any(usat) && isa(tmp, 'ADI')
                        dbo_drsu = tmp.jac{1}.diagonal(usat);
                        bou = bo(usat);
                        rsu = rs(usat);
                    else
                        doUsat = false;
                    end
                end
            else
                doUsat = false;
                rs = 0;
            end
            
            alpha = 1./(1 - dis*vap*rs.*rv);

            nph = numel(b);
            w = cell(1, nph);
            names = model.getComponentNames();
            if isempty(names)
                % No components defined?
                names = {'water', 'oil', 'gas'};
                names = names(model.getActivePhases());
            end
            pvtreg = model.PVTPropertyFunctions.Density.regions;
            for ph = 1:nph
                switch names{ph}
                    case 'water'
                        f = 1./rho{ph};
                    case 'oil'
                        if dis
                            rhoOS = rhoS(pvtreg, ph);
                            f = (alpha./rhoOS).*(1./b{oix} - dis.*rs./b{gix});
                            if doUsat
                                f(usat) = (1./(rhoOS(usat).*bou)).*(1 + (rsu./bou).*dbo_drsu);
                            end
                        else
                            f = 1./rho{ph};
                        end
                    case 'gas'
                        rhoGS = rhoS(pvtreg, ph);
                        if vap
                            f = (alpha./rhoGS).*(1./b{ph} - vap.*rv./b{oix});
                        else
                            f = 1./rho{ph};
                            if doUsat
                                f(usat) = -(1./(rhoGS(usat).*bou.^2)).*dbo_drsu;
                            end
                        end
                    otherwise
                        f = 0;
                end
                w{ph} = f./pv; % Scale final pressure equation with pore-volume
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
