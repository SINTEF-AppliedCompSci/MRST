classdef OilComponent < ImmiscibleComponent
    properties
    end
    
    methods
        function c = OilComponent(name, gasIndex)
            c@ImmiscibleComponent(name, gasIndex);
        end
        
        function c = getComponentDensity(component, model, state)
            c = getComponentDensity@ComponentImplementation(component, model, state);
            phasenames = model.getPhaseNames();
            gix = phasenames == 'G';
            oix = phasenames == 'O';
            if model.disgas
                b = model.getProps(state, 'ShrinkageFactors');
                rhoS = model.getSurfaceDensities();
                rhoOS = rhoS(oix);
                bO = b{oix};
                c{oix} = rhoOS.*bO;
            end
            assert(~model.vapoil);
        end
    end
end