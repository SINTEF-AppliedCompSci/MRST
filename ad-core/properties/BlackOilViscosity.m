classdef BlackOilViscosity < GridProperty
    properties
        useSaturatedFlag = true;
    end
    
    methods
        function mu = evaluateOnDomain(prop, model, state)
            [act, phInd] = model.getActivePhases();
            nph = sum(act);
            mu = cell(1, nph);
            
            f = model.fluid;
            [p, p_phase] = model.getProps(state, 'pressure', 'phasepressures');
            nc = numelValue(p);
            if model.water
                wix = phInd == 1;
                pw = p_phase{wix};
                mu{wix} = prop.evaluateFunctionOnGrid(f.muW, pw);
            end
            
            if model.oil
                oix = phInd == 2;
                po = p_phase{oix};
                if model.disgas
                    rs = model.getProp(state, 'rs');
                    if prop.useSaturatedFlag
                        sG = model.getProp(state, 'sg');
                        flag = sG > 0;
                    else
                        flag = false(nc, 1);
                    end
                    mu{oix} = prop.evaluateFunctionOnGrid(f.muO, po, rs, flag);
                else
                    mu{oix} = prop.evaluateFunctionOnGrid(f.muO, po);
                end
            end
            
            if model.gas
                gix = phInd == 3;
                pg = p_phase{gix};
                if model.vapoil
                    rv = model.getProp(state, 'rv');
                    if prop.useSaturatedFlag
                        sO = model.getProp(state, 'so');
                        flag = sO > 0;
                    else
                        flag = false(nc, 1);
                    end
                    mu{gix} = prop.evaluateFunctionOnGrid(f.muG, pg, rv, flag);
                else
                    mu{gix} = prop.evaluateFunctionOnGrid(f.muG, pg);
                end
            end
        end
    end
end