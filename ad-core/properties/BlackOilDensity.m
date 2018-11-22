classdef BlackOilDensity < GridProperty
    properties
    end
    
    methods
        function rho = evaluateOnGrid(prop, model, state)
            [act, phInd] = model.getActivePhases();
            fp = state.FlowProps;
            nph = sum(act);
            rho = cell(1, nph);
            
            f = model.fluid;
            p = model.getProp(state, 'pressure');
            p_phase = getPhasePressures(p, fp.CapillaryPressure);
            
            if model.water
                wix = phInd == 1;
                pw = p_phase{wix};
                bW = prop.evaluateFunctionOnGrid(f.bW, pw);
                rho{wix} = f.rhoWS.*bW;
            end
            
            if model.oil
                oix = phInd == 2;
                po = p_phase{oix};
                if model.disgas
                    rs = model.getProp(state, 'rs');
                    flag = false(size(double(p)));
                    bO = prop.evaluateFunctionOnGrid(f.bO, po, rs, flag);
                else
                    bO = prop.evaluateFunctionOnGrid(f.bO, po);
                end
                rho{oix} = f.rhoOS.*bO;
            end
            
            if model.gas
                gix = phInd == 3;
                pg = p_phase{gix};
                if model.vapoil
                    rv = model.getProp(state, 'rv');
                    flag = false(size(double(p)));
                    bG = prop.evaluateFunctionOnGrid(f.bG, pg, rv, flag);
                else
                    bG = prop.evaluateFunctionOnGrid(f.bG, pg);
                end
                rho{gix} = f.rhoGS.*bG;
            end
        end
    end
end