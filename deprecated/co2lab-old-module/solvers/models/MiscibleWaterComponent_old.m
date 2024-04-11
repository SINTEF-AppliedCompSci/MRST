classdef MiscibleWaterComponent < ImmiscibleComponent
    
    methods
        function c = MiscibleWaterComponent(name, waterIndex)
            c@ImmiscibleComponent(name, waterIndex);
            c = c.functionDependsOn('getComponentDensity', ...
                                    {'ShrinkageFactors', 'SurfaceDensity'}, ...
                                    'PVTPropertyFunctions');
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

                wix = (phasenames == 'W');
                pvt = model.PVTPropertyFunctions;
                
                rhoS = pvt.get(model, state, 'SurfaceDensity', true);
                b = pvt.get(model, state, 'ShrinkageFactors', true);
                
                c{wix} = rhoS{wix} .* b{wix};

                
            else
                % density of water in water phase (which is the only phase
                % with water)
                c = getComponentDensity@ImmiscibleComponent(component, model, state); 
            end
        end
            
    end

    
    
end
