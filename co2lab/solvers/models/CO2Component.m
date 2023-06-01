classdef CO2Component < ImmiscibleComponent

    methods
        function c = CO2Component(name, disgas, gasIndex)
            c@ImmiscibleComponent(name, gasIndex);
            c = c.functionDependsOn('getComponentDensity', ...
                                    {'ShrinkageFactors', 'SurfaceDensity'}, ...
                                    'PVTPropertyFunctions');
            if disgas
                c = c.functionDependsOn('getComponentDensity', 'rs', 'state');
            end
        end
        
        function c = getPhaseComposition(component, model, state, varargin)
            % @@ This is modeled after the GasComponent in the black oil module.  But should 
            %    it not rather call the GenericComponent version, to account for presence in 
            %    multiple phases?
            c = getPhaseComposition@ImmiscibleComponent(component, model, state, varargin{:});
        end
        
        function c = getComponentDensity(component, model, state, varargin)
            if model.disgas
                % establish (empty) cell array with one entry per phase
                phasenames = model.getPhaseNames();
                nph = numel(phasenames);
                c = cell(nph, 1);
                
                gix = (phasenames == 'G');
                pvt = model.PVTPropertyFunctions;
                
                rho = pvt.get(model, state, 'Density', true);
                b = pvt.get(model, state, 'ShrinkageFactors', true);
                rhoS = pvt.get(model, state, 'SurfaceDensity', true);
                rhoGS = rhoS{gix};
                
                % density of CO2 in CO2 phase
                c{gix} = rho{gix};
                
                % density of CO2 in brine phase
                if model.disgas
                    wix = (phasenames == 'W');
                    bW = b{wix};
                    rs = model.getProp(state, 'rs');
                    c{wix} = rs .* rhoGS .* bW;
                end
            else
                % density of CO2 in CO2 phase (which is the only phase with CO2)
                c = getComponentDensity@ImmiscibleComponent(component, model, state, varargin{:});
            end
        end
    end
end
