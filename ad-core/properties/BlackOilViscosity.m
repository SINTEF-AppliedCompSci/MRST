classdef BlackOilViscosity < GridProperty
    properties
    end
    
    methods
        function mu = evaluateOnGrid(prop, model, state)
            [act, phInd] = model.getActivePhases();
            fp = state.FlowProps;
            nph = sum(act);
            mu = cell(1, nph);
            
            f = model.fluid;
            p = model.getProp(state, 'pressure');
            if model.water
                wix = phInd == 1;
                pw = p;
                pcwo = fp.CapillaryPressure{wix};
                if ~isempty(pcwo)
                    pw = pw + pcwo;
                end
                mu{wix} = prop.evaluateFunctionOnGrid(f.muW, pw);
            end
            
            if model.oil
                oix = phInd == 2;
                if model.disgas
                    rs = model.getProp(state, 'rs');
                    flag = false(size(double(p)));
                    mu{oix} = prop.evaluateFunctionOnGrid(f.muO, p, rs, flag);
                else
                    mu{oix} = prop.evaluateFunctionOnGrid(f.muO, p);
                end
            end
            
            if model.gas
                gix = phInd == 3;
                pg = p;
                pcgo = fp.CapillaryPressure{gix};
                if ~isempty(pcwo)
                    pg = pg + pcgo;
                end
                if model.vapoil
                    rv = model.getProp(state, 'rv');
                    flag = false(size(double(p)));
                    mu{gix} = prop.evaluateFunctionOnGrid(f.muG, pg, rv, flag);
                else
                    mu{gix} = prop.evaluateFunctionOnGrid(f.muG, pg);
                end
            end
        end
    end
end