classdef BlackOilViscosity < StateFunction
    properties
        useSaturatedFlag = false;
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
        end
        
        function mu = evaluateOnDomain(prop, model, state)
            [act, phInd] = model.getActivePhases();
            nph = sum(act);
            mu = cell(1, nph);
            
            f = model.fluid;
            p_phase = prop.getEvaluatedDependencies(state, 'PhasePressures');
            nc = numelValue(p_phase{1});
            if model.water
                wix = phInd == 1;
                pw = p_phase{wix};
                mu{wix} = prop.evaluateFunctionOnDomainWithArguments(f.muW, pw);
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
                    mu{oix} = prop.evaluateFunctionOnDomainWithArguments(f.muO, po, rs, flag);
                else
                    mu{oix} = prop.evaluateFunctionOnDomainWithArguments(f.muO, po);
                end
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
                    mu{gix} = prop.evaluateFunctionOnDomainWithArguments(f.muG, pg, rv, flag);
                else
                    mu{gix} = prop.evaluateFunctionOnDomainWithArguments(f.muG, pg);
                end
            end
        end
    end
end