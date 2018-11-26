classdef BlackOilDensity < GridProperty
    properties
    end
    
    methods
        function rho = evaluateOnGrid(prop, model, state)
            rhoS = model.getSurfaceDensities();
            fp = state.FlowProps;
            nph = numel(rhoS);
            rho = cell(1, nph);
            b = fp.ShrinkageFactors;
            
            for i = 1:numel(b)
                rho{i} = rhoS(i).*b{i};
            end
            if (model.disgas || model.vapoil) && model.gas && model.oil
                names = model.getPhaseNames();
                oix = names == 'O';
                gix = names == 'G';
                
                if model.oil && model.disgas
                    rs = model.getProp(state, 'rs');
                    rho{oix} = rho{oix} + rs.*b{oix}.*rhoS(gix);
                end
                if model.oil && model.disgas
                    rv = model.getProp(state, 'rv');
                    rho{gix} = rho{gix} + rv.*b{gix}.*rhoS(oix);
                end
            end
        end
    end
end