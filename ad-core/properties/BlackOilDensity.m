classdef BlackOilDensity < GridProperty
    properties
    end
    
    methods
        function rho = evaluateOnDomain(prop, model, state)
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
                
                if model.disgas
                    rs = model.getProp(state, 'rs');
                    rho{oix} = rho{oix} + rs.*b{oix}.*rhoS(gix);
                end
                if model.vapoil
                    rv = model.getProp(state, 'rv');
                    rho{gix} = rho{gix} + rv.*b{gix}.*rhoS(oix);
                end
            end
            
            mv = cellfun(@(x) min(value(x)), rho);
            if any(mv <= 0)
            	warning('Negative densities detected! Capping to 1e-12.')
                rho = cellfun(@(x) max(x, 1e-12), rho, 'UniformOutput', false);
            end
        end
    end
end