classdef BlackOilViscosity < GridProperty
    properties
    end
    
    methods
        function rho = evaluateOnGrid(prop, model, state)
            act = model.getActivePhases();
            nph = sum(act);
            rho = cell(1, nph);
            ix = 1;
            
            f = model.fluid;
            p = model.getProp(state, 'pressure');
            if model.water
                bW = prop.evaluateFunctionOnGrid(f.bW, p);
                rho{ix} = f.rhoWS.*bW;
                ix = ix + 1;
            end
            
            if model.oil
                if model.disgas
                    rs = model.getProp(state, 'rs');
                    flag = false(size(double(p)));
                    bO = prop.evaluateFunctionOnGrid(f.bO, p, rs, flag);
                else
                    bO = prop.evaluateFunctionOnGrid(f.bO, p);
                end
                rho{ix} = f.rhoOS.*bO;
                ix = ix + 1;
            end
            
            if model.gas
                if model.vapoil
                    rv = model.getProp(state, 'rv');
                    flag = false(size(double(p)));
                    bG = prop.evaluateFunctionOnGrid(f.bG, p, rv, flag);
                else
                    bG = prop.evaluateFunctionOnGrid(f.bG, p);
                end
                rho{ix} = f.rhoGS.*bG;
            end
        end
    end
end