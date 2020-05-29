classdef BlackOilPressureReductionFactors < PressureReductionFactors
    % Component weighting factors used to form a pressure equation
    properties
        disgas = false;
        vapoil = false;
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
            
            [act, phInd] = model.getActivePhases();
            if vap
                rv = model.getProp(state, 'rv');
                oix = phInd == 2;
            else
                rv = 0;
            end
            
            if dis
                rs = model.getProp(state, 'rs');
                gix = phInd == 3;
            else
                rs = 0;
            end
            
            alpha = 1./(1 - dis*vap*rs.*rv);

            nph = sum(act);
            w = cell(1, nph);
            names = model.getComponentNames();
            for ph = 1:nph
                switch names{ph}
                    case 'water'
                        f = 1./rho{ph};
                    case 'oil'
                        if dis
                            f = (alpha./rhoS(ph)).*(1./b{ph} - dis.*rs./b{gix});
                        else
                            f = 1./rho{ph};
                        end
                    case 'gas'
                        if vap
                            f = (alpha./rhoS(ph)).*(1./b{ph} - vap.*rv./b{oix});
                        else
                            f = 1./rho{ph};
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
Copyright 2009-2020 SINTEF Digital, Mathematics & Cybernetics.

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
