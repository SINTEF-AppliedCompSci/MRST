classdef ConnectionPressureDrop < StateFunction
   
    properties
        
    end
    
    methods
        
        function gp = ConnectionPressureDrop(varargin)
            
            gp@StateFunction(varargin{:});
            gp = gp.dependsOn({'FacilityWellMapping'});
            gp = gp.dependsOn('Density', 'PVTPropertyFunctions');
            gp = gp.dependsOn('bhp', 'state');
            gp.label = 'p_c-p_{bh}-g \Delta z \rho_{w}';
            
        end
        
        function cdp = evaluateOnDomain(prop, model, state)
            
            
            map = prop.getEvaluatedDependencies(state, 'FacilityWellMapping');
           
            [rho] = model.ReservoirModel.getProps(state, 'Density');
            rho = cellfun(@(rho) rho(map.cells), rho, 'UniformOutput', false);
            rho = rho{1};
            g = norm(model.ReservoirModel.gravity);
            
            nw = numel(map.W);
            cdp = cell(nw,1);
            
            for i = 1:nw
                dz = diff([0; map.W(i).dZ]);
                rhoW = rho(map.perf2well == i);
                np = numel(dz);
                S = ones(np);
                S = sparse(tril(S));
                cdpW = S*(g.*rhoW.*dz);
                cdp{i} = cdpW;
            end
            cdp = vertcat(cdp{:});
            
        end
        
    end
    
end