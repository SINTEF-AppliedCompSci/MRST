classdef BlackOilDensity < StateFunction
    properties
        disgas = false;
        vapoil = false;
    end
    
    methods
        function gp = BlackOilDensity(model, varargin)
            gp@StateFunction(model, varargin{:});
            if isprop(model, 'disgas')
                gp.disgas = model.disgas;
                gp = gp.dependsOn({'rs'}, 'state');
            end
            if isprop(model, 'vapoil')
                gp.vapoil = model.vapoil;
                gp = gp.dependsOn({'rv'}, 'state');
            end
            gp = gp.dependsOn({'ShrinkageFactors'});
        end
        function rho = evaluateOnDomain(prop, model, state)
            rhoS = model.getSurfaceDensities();
            nph = size(rhoS, 2);
            rho = cell(1, nph);
            b = prop.getEvaluatedDependencies(state, 'ShrinkageFactors');
            for i = 1:numel(b)
                rho{i} = rhoS(i).*b{i};
            end
            if (prop.disgas || prop.vapoil) && model.gas && model.oil
                names = model.getPhaseNames();
                oix = names == 'O';
                gix = names == 'G';
                
                if prop.disgas
                    rs = model.getProp(state, 'rs');
                    rho{oix} = rho{oix} + rs.*b{oix}.*rhoS(prop.regions, gix);
                end
                if prop.vapoil
                    rv = model.getProp(state, 'rv');
                    rho{gix} = rho{gix} + rv.*b{gix}.*rhoS(prop.regions, oix);
                end
            end
            
            mv = cellfun(@(x) min(value(x)), rho);
            if model.verbose > 1 && any(mv <= 0)
            	warning('Negative densities detected! Capping to 1e-12.')
                rho = cellfun(@(x) max(x, 1e-12), rho, 'UniformOutput', false);
            end
        end
    end
end