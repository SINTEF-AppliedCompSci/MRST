classdef BlackOilShrinkageFactors < GridProperty
    properties
        useSaturatedFlag = false;
        disgas = false;
        vapoil = false;
    end
    
    methods
        function gp = BlackOilShrinkageFactors(model, varargin)
            gp@GridProperty(model, varargin{:});
            if isprop(model, 'disgas')
                gp.disgas = model.disgas;
            end
            if isprop(model, 'vapoil')
                gp.vapoil = model.vapoil;
            end
        end

        function b = evaluateOnDomain(prop, model, state)
            [act, phInd] = model.getActivePhases();
            nph = sum(act);
            b = cell(1, nph);
            
            f = model.fluid;
            [p, p_phase] = model.getProps(state, 'pressure', 'phasepressures');
            
            
            if model.water
                wix = phInd == 1;
                pw = p_phase{wix};
                bW = prop.evaluateFunctionOnGrid(f.bW, pw);
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
                        flag = false(numelValue(p), 1);
                    end
                    bO = prop.evaluateFunctionOnGrid(f.bO, po, rs, flag);
                else
                    bO = prop.evaluateFunctionOnGrid(f.bO, po);
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
                        flag = false(numelValue(p), 1);
                    end
                    bG = prop.evaluateFunctionOnGrid(f.bG, pg, rv, flag);
                else
                    bG = prop.evaluateFunctionOnGrid(f.bG, pg);
                end
                b{gix} = bG;
            end
        end
    end
end